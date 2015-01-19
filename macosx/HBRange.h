/*  HBRange.h $

 This file is part of the HandBrake source code.
 Homepage: <http://handbrake.fr/>.
 It may be used under the terms of the GNU General Public License. */

#import <Foundation/Foundation.h>

@class HBTitle;

typedef NS_ENUM(NSUInteger, HBRangeType) {
    HBRangeTypeChapters,
    HBRangeTypeFrames,
    HBRangeTypeSeconds,
    HBRangePreviewIndex,
};

@interface HBRange : NSObject <NSCoding, NSCopying>

- (instancetype)initWithTitle:(HBTitle *)title;

@property (nonatomic, readwrite) HBRangeType type;

/// HBRangeTypeChapters
@property (nonatomic, readwrite) int chapterStart;
@property (nonatomic, readwrite) int chapterStop;

/// HBRangeTypeFrames
@property (nonatomic, readwrite) int frameStart;
@property (nonatomic, readwrite) int frameStop;

/// HBRangeTypeSeconds
@property (nonatomic, readwrite) int secondsStart;
@property (nonatomic, readwrite) int secondsStop;

/// HBRangePreviewIndex
@property (nonatomic, readwrite) int previewIndex;
@property (nonatomic, readwrite) int previewsCount;
@property (nonatomic, readwrite) int64_t ptsToStop;

@property (nonatomic, readonly) NSString *duration;

@property (nonatomic, readwrite, assign) HBTitle *title;

@end
