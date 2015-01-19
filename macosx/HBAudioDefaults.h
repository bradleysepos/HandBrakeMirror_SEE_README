/*  HBAudioSettings.h $

 This file is part of the HandBrake source code.
 Homepage: <http://handbrake.fr/>.
 It may be used under the terms of the GNU General Public License. */

#import <Foundation/Foundation.h>
#import "HBPresetCoding.h"

typedef NS_ENUM(NSUInteger, HBAudioTrackSelectionBehavior) {
    HBAudioTrackSelectionBehaviorNone,
    HBAudioTrackSelectionBehaviorFirst,
    HBAudioTrackSelectionBehaviorAll,
};

/**
 *  HBAudioSettings
 *  Stores the audio defaults settings.
 */
@interface HBAudioDefaults : NSObject <NSCoding, NSCopying, HBPresetCoding>

@property (nonatomic, readwrite) HBAudioTrackSelectionBehavior trackSelectionBehavior;
@property (nonatomic, readwrite, retain) NSMutableArray *trackSelectionLanguages;

@property (nonatomic, readwrite, retain) NSMutableArray *tracksArray;

@property(nonatomic, readwrite) BOOL allowAACPassthru;
@property(nonatomic, readwrite) BOOL allowAC3Passthru;
@property(nonatomic, readwrite) BOOL allowDTSHDPassthru;
@property(nonatomic, readwrite) BOOL allowDTSPassthru;
@property(nonatomic, readwrite) BOOL allowMP3Passthru;

@property(nonatomic, readwrite) int encoderFallback;
@property(nonatomic, readwrite) BOOL secondaryEncoderMode;

@property(nonatomic, readonly) NSArray *audioEncoderFallbacks;

- (void)validateEncoderFallbackForVideoContainer:(int)container;

@end
