/*  HBSubtitlesSettings.m $

 This file is part of the HandBrake source code.
 Homepage: <http://handbrake.fr/>.
 It may be used under the terms of the GNU General Public License. */

#import "HBSubtitlesDefaults.h"
#import "NSCodingMacro.h"

@implementation HBSubtitlesDefaults

- (instancetype)init
{
    self = [super init];
    if (self)
    {
        _trackSelectionLanguages = [[NSMutableArray alloc] init];
    }
    return self;
}

- (void)applyPreset:(NSDictionary *)preset
{
    if ([preset[@"SubtitleTrackSelectionBehavior"] isEqualToString:@"first"])
    {
        self.trackSelectionBehavior = HBSubtitleTrackSelectionBehaviorFirst;
    }
    else if ([preset[@"SubtitleTrackSelectionBehavior"] isEqualToString:@"all"])
    {
        self.trackSelectionBehavior = HBSubtitleTrackSelectionBehaviorAll;
    }
    else
    {
        self.trackSelectionBehavior = HBSubtitleTrackSelectionBehaviorNone;
    }
    self.trackSelectionLanguages = [NSMutableArray arrayWithArray:preset[@"SubtitleLanguageList"]];
    self.addCC = [preset[@"SubtitleAddCC"] boolValue];
    self.addForeignAudioSearch = [preset[@"SubtitleAddForeignAudioSearch"] boolValue];
    self.addForeignAudioSubtitle = [preset[@"SubtitleAddForeignAudioSubtitle"] boolValue];
}

- (void)writeToPreset:(NSMutableDictionary *)preset
{
    if (self.trackSelectionBehavior == HBSubtitleTrackSelectionBehaviorFirst)
    {
        preset[@"SubtitleTrackSelectionBehavior"] = @"first";
    }
    else if (self.trackSelectionBehavior == HBSubtitleTrackSelectionBehaviorAll)
    {
        preset[@"SubtitleTrackSelectionBehavior"] = @"all";
    }
    else
    {
        preset[@"SubtitleTrackSelectionBehavior"] = @"none";
    }

    preset[@"SubtitleLanguageList"] = [[self.trackSelectionLanguages copy] autorelease];
    preset[@"SubtitleAddCC"] = @(self.addCC);
    preset[@"SubtitleAddForeignAudioSearch"] = @(self.addForeignAudioSearch);
    preset[@"SubtitleAddForeignAudioSubtitle"] = @(self.addForeignAudioSubtitle);
}

#pragma mark - NSCopying

- (instancetype)copyWithZone:(NSZone *)zone
{
    HBSubtitlesDefaults *copy = [[[self class] alloc] init];

    if (copy)
    {
        copy->_trackSelectionBehavior = _trackSelectionBehavior;
        [copy->_trackSelectionLanguages release];
        copy->_trackSelectionLanguages = [_trackSelectionLanguages mutableCopy];

        copy->_addForeignAudioSearch = _addForeignAudioSearch;
        copy->_addForeignAudioSubtitle = _addForeignAudioSubtitle;
        copy->_addCC = _addCC;
    }
    
    return copy;
}

#pragma mark - NSCoding

- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeInt:1 forKey:@"HBAudioDefaultsVersion"];

    encodeInteger(_trackSelectionBehavior);
    encodeObject(_trackSelectionLanguages);

    encodeBool(_addForeignAudioSearch);
    encodeBool(_addForeignAudioSubtitle);
    encodeBool(_addCC);
}

- (id)initWithCoder:(NSCoder *)decoder
{
    self = [super init];

    decodeInteger(_trackSelectionBehavior);
    decodeObject(_trackSelectionLanguages);

    decodeBool(_addForeignAudioSearch);
    decodeBool(_addForeignAudioSubtitle);
    decodeBool(_addCC);

    return self;
}

- (void)dealloc
{
    [_trackSelectionLanguages release];
    [super dealloc];
}

@end
