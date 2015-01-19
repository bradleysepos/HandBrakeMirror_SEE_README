#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <glib.h>
#include <glib/gstdio.h>
#include "plist.h"
#include "values.h"
#include <gdk-pixbuf/gdk-pixbuf.h>
#include <gdk-pixbuf/gdk-pixdata.h>
#include <gtk/gtk.h>

enum
{
    R_NONE = 0,
    R_RESOURCES,
    R_SECTION,
    R_ICON,
    R_PLIST,
    R_STRING,
};

typedef struct
{
    gchar *tag;
    gint id;
} tag_map_t;

static tag_map_t tag_map[] =
{
    {"resources", R_RESOURCES},
    {"section", R_SECTION},
    {"icon", R_ICON},
    {"plist", R_PLIST},
    {"string", R_STRING},
};
#define TAG_MAP_SZ  (sizeof(tag_map)/sizeof(tag_map_t))

typedef struct
{
    gchar *key;
    gchar *value;
    GValue *plist;
    GQueue *stack;
    GQueue *tag_stack;
    gboolean closed_top;
} parse_data_t;

GList *inc_list = NULL;

static gchar*
find_file(GList *list, const gchar *name)
{
    gchar *str;
    GList *link = list;

    while (link != NULL)
    {
        gchar *inc;

        inc = (gchar*)link->data;
        str = g_strdup_printf("%s/%s", inc, name);
        if (g_file_test(str, G_FILE_TEST_IS_REGULAR))
        {
            return str;
        }
        g_free(str);
        link = g_list_next(link);
    }
    if (g_file_test(name, G_FILE_TEST_IS_REGULAR))
    {
        return g_strdup(name);
    }
    return NULL;
}

static const gchar*
lookup_attr_value(
    const gchar *name,
    const gchar **attr_names,
    const gchar **attr_values)
{
    gint ii;

    if (attr_names == NULL) return NULL;
    for (ii = 0; attr_names[ii] != NULL; ii++)
    {
        if (strcmp(name, attr_names[ii]) == 0)
            return attr_values[ii];
    }
    return NULL;
}

static GValue*
read_string_from_file(const gchar *fname)
{
    gchar *buffer;
    size_t size;
    GValue *gval;
    FILE *fd;

    fd = g_fopen(fname, "r");
    if (fd == NULL)
        return NULL;
    fseek(fd, 0, SEEK_END);
    size = ftell(fd);
    fseek(fd, 0, SEEK_SET);
    buffer = g_malloc(size+1);
    size = fread(buffer, 1, size, fd);
    buffer[size] = 0;
    gval = ghb_value_new(G_TYPE_STRING);
    g_value_take_string(gval, buffer);
    fclose(fd);
    return gval;
}

static void add_icon(GValue *dict, const char *fname)
{
    FILE *f;

    GdkPixbufFormat *pbf;
    int width, height;
    gboolean svg;

    pbf = gdk_pixbuf_get_file_info(fname, &width, &height);
    svg = gdk_pixbuf_format_is_scalable(pbf);

    f = fopen(fname, "rb");
    if (f == NULL)
    {
        fprintf(stderr, "open failed: %s\n", fname);
        return;
    }

    ghb_rawdata_t *rd;
    rd = g_malloc(sizeof(ghb_rawdata_t));
    fseek(f, 0, SEEK_END);
    rd->size = ftell(f);
    fseek(f, 0, SEEK_SET);

    rd->data = g_malloc(rd->size);
    fread(rd->data, 1, rd->size, f);
    
    GValue *data = ghb_rawdata_value_new(rd);
    ghb_dict_insert(dict, g_strdup("svg"), ghb_boolean_value_new(svg));
    ghb_dict_insert(dict, g_strdup("data"), data);
}

static void insert_value(GValue *container, const char *key, GValue *element)
{
    GType gtype;

    gtype = G_VALUE_TYPE(container);
    if (gtype == ghb_array_get_type())
    {
        ghb_array_append(container, element);
    }
    else if (gtype == ghb_dict_get_type())
    {
        if (key == NULL)
        {
            g_warning("No key for dictionary item");
            ghb_value_free(element);
        }
        else
        {
            ghb_dict_insert(container, g_strdup(key), element);
        }
    }
    else
    {
        g_error("Invalid container type. This shouldn't happen");
    }
}

static void
start_element(
    GMarkupParseContext *ctx,
    const gchar *tag,
    const gchar **attr_names,
    const gchar **attr_values,
    gpointer ud,
    GError **error)
{
    parse_data_t *pd = (parse_data_t*)ud;
    union
    {
        gint id;
        gpointer pid;
    } id;
    gint ii;

    // Check to see if the first element found has been closed
    // If so, ignore any junk following it.
    if (pd->closed_top)
        return;

    for (ii = 0; ii < TAG_MAP_SZ; ii++)
    {
        if (strcmp(tag, tag_map[ii].tag) == 0)
        {
            id.id = tag_map[ii].id;
            break;
        }
    }
    if (ii == TAG_MAP_SZ)
    {
        g_warning("Unrecognized start tag (%s)", tag);
        return;
    }
    g_queue_push_head(pd->tag_stack, id.pid);
    GType gtype = 0;
    GValue *gval = NULL;
    GValue *current = g_queue_peek_head(pd->stack);
    switch (id.id)
    {
        case R_SECTION:
        {
            const gchar *key;

            key = lookup_attr_value("name", attr_names, attr_values);
            if (key == NULL)
            {
                g_warning("section: missing a requried *name* attribute");
                exit(EXIT_FAILURE);
            }
            if (strcmp(key, "icons") == 0)
            {
                gval = ghb_dict_value_new();
                if (pd->key) g_free(pd->key);
                pd->key = g_strdup(key);
                g_queue_push_head(pd->stack, gval);
            }
        } break;
        case R_ICON:
        {
            gchar *fname;
            const gchar *name, *key;

            name = lookup_attr_value("file", attr_names, attr_values);
            if (name == NULL)
            {
                g_warning("icon: missing a requried *file* attribute");
                exit(EXIT_FAILURE);
            }
            fname = find_file(inc_list, name);
            if (fname == NULL)
            {
                g_warning("icon: no such file %s", name);
                exit(EXIT_FAILURE);
            }
            key = lookup_attr_value("name", attr_names, attr_values);
            if (key == NULL)
            {
                g_warning("icon: missing a requried *name* attribute");
                g_free(fname);
                exit(EXIT_FAILURE);
            }
            gval = ghb_dict_value_new();
            add_icon(gval, fname);

            g_free(pd->key);
            pd->key = g_strdup(key);
            g_free(fname);
        } break;
        case R_PLIST:
        {
            gchar *fname;
            const gchar *name, *key;

            name = lookup_attr_value("file", attr_names, attr_values);
            if (name == NULL)
            {
                g_warning("plist: missing a requried *file* attribute");
                exit(EXIT_FAILURE);
            }
            fname = find_file(inc_list, name);
            if (fname == NULL)
            {
                g_warning("plist: no such file %s", name);
                exit(EXIT_FAILURE);
            }
            key = lookup_attr_value("name", attr_names, attr_values);
            if (key == NULL)
            {
                g_warning("plist: missing a requried *name* attribute");
                g_free(fname);
                exit(EXIT_FAILURE);
            }
            gval = ghb_plist_parse_file(fname);
            if (pd->key) g_free(pd->key);
            pd->key = g_strdup(key);
            g_free(fname);
        } break;
        case R_STRING:
        {
            gchar *fname;
            const gchar *name, *key;
            const gchar *version;
            char *end;
            int major = 0, minor = 0, micro = 0;

            name = lookup_attr_value("file", attr_names, attr_values);
            if (name == NULL)
            {
                g_warning("string: missing a requried *file* attribute");
                exit(EXIT_FAILURE);
            }
            fname = find_file(inc_list, name);
            if (fname == NULL)
            {
                g_warning("string: no such file %s", name);
                exit(EXIT_FAILURE);
            }
            key = lookup_attr_value("name", attr_names, attr_values);
            if (key == NULL)
            {
                g_warning("string: missing a requried *name* attribute");
                g_free(fname);
                exit(EXIT_FAILURE);
            }
            version = lookup_attr_value("version", attr_names, attr_values);
            if (version != NULL)
            {
                major = strtol(version, &end, 10);
                if (end != version && *end != 0)
                {
                    version = end + 1;
                    minor = strtol(version, &end, 10);
                    if (end != version && *end != 0)
                    {
                        version = end + 1;
                        micro = strtol(version, &end, 10);
                        if (end != version && *end != 0)
                        {
                            version = end + 1;
                        }
                    }
                }
            }
            if (GTK_CHECK_VERSION(major, minor, micro))
            {
                gval = read_string_from_file(fname);
                if (pd->key) g_free(pd->key);
                pd->key = g_strdup(key);
            }
            g_free(fname);
        } break;
    }
    // Add the element to the current container
    if (gval)
    { // There's an element to add
        if (current == NULL)
        {
            pd->plist = gval;
            return;
        }
        insert_value(current, pd->key, gval);
    }
}

static void
end_element(
    GMarkupParseContext *ctx,
    const gchar *name,
    gpointer ud,
    GError **error)
{
    parse_data_t *pd = (parse_data_t*)ud;
    gint id;
    union
    {
        gint id;
        gpointer pid;
    } start_id;
    gint ii;

    // Check to see if the first element found has been closed
    // If so, ignore any junk following it.
    if (pd->closed_top)
        return;

    for (ii = 0; ii < TAG_MAP_SZ; ii++)
    {
        if (strcmp(name, tag_map[ii].tag) == 0)
        {
            id = tag_map[ii].id;
            break;
        }
    }
    if (ii == TAG_MAP_SZ)
    {
        g_warning("Unrecognized start tag (%s)", name);
        return;
    }
    start_id.pid = g_queue_pop_head(pd->tag_stack);
    if (start_id.id != id)
        g_warning("start tag != end tag: (%s %d) %d", name, id, id);

    GValue *gval = NULL;
    GValue *current = g_queue_peek_head(pd->stack);
    GType gtype = 0;
    switch (id)
    {
        case R_SECTION:
        {
            g_queue_pop_head(pd->stack);
        } break;
    }
    if (gval)
    {
        // Get the top of the data structure stack and if it's an array
        // or dict, add the current element
        if (current == NULL)
        {
            pd->plist = gval;
            pd->closed_top = TRUE;
            return;
        }
        insert_value(current, pd->key, gval);
    }
    if (g_queue_is_empty(pd->tag_stack))
        pd->closed_top = TRUE;
}

static void
text_data(
    GMarkupParseContext *ctx,
    const gchar *text,
    gsize len,
    gpointer ud,
    GError **error)
{
    parse_data_t *pd = (parse_data_t*)ud;
    if (pd->value) g_free(pd->value);
    pd->value = g_strdup(text);
}

static void
passthrough(
    GMarkupParseContext *ctx,
    const gchar *text,
    gsize len,
    gpointer ud,
    GError **error)
{
    //parse_data_t *pd = (parse_data_t*)ud;

    //g_debug("passthrough %s", text);
}

static void
parse_error(GMarkupParseContext *ctx, GError *error, gpointer ud)
{
    g_warning("Resource parse error: %s", error->message);
}

// This is required or the parser crashes
static void
destroy_notify(gpointer data)
{ // Do nothing
    //g_debug("destroy parser");
}

GValue*
ghb_resource_parse(const gchar *buf, gssize len)
{
    GMarkupParseContext *ctx;
    GMarkupParser parser;
    parse_data_t pd;
    GError *err = NULL;

    pd.stack = g_queue_new();
    pd.tag_stack = g_queue_new();
    pd.key = NULL;
    pd.value = NULL;
    pd.plist = ghb_dict_value_new();
    g_queue_push_head(pd.stack, pd.plist);
    pd.closed_top = FALSE;

    parser.start_element = start_element;
    parser.end_element = end_element;
    parser.text = text_data;
    parser.passthrough = passthrough;
    parser.error = parse_error;
    ctx = g_markup_parse_context_new(&parser, 0, &pd, destroy_notify);

    g_markup_parse_context_parse(ctx, buf, len, &err);
    g_markup_parse_context_end_parse(ctx, &err);
    g_markup_parse_context_free(ctx);
    g_queue_free(pd.stack);
    g_queue_free(pd.tag_stack);
    return pd.plist;
}

GValue*
ghb_resource_parse_file(FILE *fd)
{
    gchar *buffer;
    size_t size;
    GValue *gval;

    if (fd == NULL)
        return NULL;
    fseek(fd, 0, SEEK_END);
    size = ftell(fd);
    fseek(fd, 0, SEEK_SET);
    buffer = g_malloc(size+1);
    size = fread(buffer, 1, size, fd);
    buffer[size] = 0;
    gval = ghb_resource_parse(buffer, (gssize)size);
    g_free(buffer);
    return gval;
}

static void
usage(char *cmd)
{
    fprintf(stderr,
"Usage: %s [-I <inc path>] <in resource list> <out resource plist>\n"
"Summary:\n"
"    Creates a resource plist from a resource list\n"
"Options:\n"
"    I - Include path to search for files\n"
"    <in resource list>    Input resources file\n"
"    <out resource plist>  Output resources plist file\n"
, cmd);

    exit(EXIT_FAILURE);
}

#define OPTS "I:"

gint
main(gint argc, gchar *argv[])
{
    FILE *file;
    GValue *gval;
    int opt;
    const gchar *src, *dst;

    do
    {
        opt = getopt(argc, argv, OPTS);
        switch (opt)
        {
        case -1: break;

        case 'I':
            inc_list = g_list_prepend(inc_list, g_strdup(optarg));
            break;
        }
    } while (opt != -1);

    if (optind != argc - 2)
    {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    src = argv[optind++];
    dst = argv[optind++];

#if !GLIB_CHECK_VERSION(2, 36, 0)
    g_type_init();
#endif

    file = g_fopen(src, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open %s\n", src);
        return EXIT_FAILURE;
    }

    gval = ghb_resource_parse_file(file);
    ghb_plist_write_file(dst, gval);
    fclose(file);
    return 0;
}

