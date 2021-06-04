
#include <cstdint>
#include <stdlib.h> 
#include <assert.h>
#include <iostream>
using namespace std;

#define KVZ_BIT_DEPTH 8
#define SWAP(a,b,swaptype) { swaptype tempval; tempval = a; a = b; b = tempval; }
#define MAX(a,b) (((a)>(b))?(a):(b))
#define PIXEL_MAX ((1 << KVZ_BIT_DEPTH) - 1)
#define CLIP(low,high,value) MAX((low),MIN((high),(value)))
#define CLIP_TO_PIXEL(value) CLIP(0, PIXEL_MAX, (value))

typedef enum { COLOR_Y = 0, COLOR_U, COLOR_V } color_t;
typedef uint8_t kvz_pixel;
typedef struct {
    kvz_pixel left[2 * 32 + 1];
    kvz_pixel top[2 * 32 + 1];
} kvz_intra_ref;
typedef struct
{
    kvz_intra_ref ref;
    kvz_intra_ref filtered_ref;
    bool filtered_initialized;
} kvz_intra_references;


inline static void intra_filter_reference(
    int_fast8_t log2_width,
    kvz_intra_references* refs)
{
    if (refs->filtered_initialized) {
        return;
    }
    else {
        refs->filtered_initialized = true;
    }

    const int_fast8_t ref_width = 2 * (1 << log2_width) + 1;
    kvz_intra_ref* ref = &refs->ref;
    kvz_intra_ref* filtered_ref = &refs->filtered_ref;

    filtered_ref->left[0] = (ref->left[1] + 2 * ref->left[0] + ref->top[1] + 2) / 4;
    filtered_ref->top[0] = filtered_ref->left[0];

    for (int_fast8_t y = 1; y < ref_width - 1; ++y) {
        kvz_pixel* p = &ref->left[y];
        filtered_ref->left[y] = (p[-1] + 2 * p[0] + p[1] + 2) / 4;
    }
    filtered_ref->left[ref_width - 1] = ref->left[ref_width - 1];

    for (int_fast8_t x = 1; x < ref_width - 1; ++x) {
        kvz_pixel* p = &ref->top[x];
        filtered_ref->top[x] = (p[-1] + 2 * p[0] + p[1] + 2) / 4;
    }
    filtered_ref->top[ref_width - 1] = ref->top[ref_width - 1];
}

inline unsigned kvz_math_floor_log2(unsigned value)
{
    assert(value > 0);

    unsigned result = 0;

    for (int i = 4; i >= 0; --i) {
        unsigned bits = 1ull << i;
        unsigned shift = value >= (1 << bits) ? bits : 0;
        result += shift;
        value >>= shift;
    }

    return result;
}
#define MIN(a,b) (((a)<(b))?(a):(b))

inline void intra_pred_dc(
    const int_fast8_t log2_width,
    const kvz_pixel* const ref_top,
    const kvz_pixel* const ref_left,
    kvz_pixel* const out_block)
{
    int_fast8_t width = 1 << log2_width;

    int_fast16_t sum = 0;
    for (int_fast8_t i = 0; i < width; ++i) {
        sum += ref_top[i + 1];
        sum += ref_left[i + 1];
    }

    const kvz_pixel dc_val = (sum + width) >> (log2_width + 1);
    const int_fast16_t block_size = 1 << (log2_width * 2);
    cout << "dc_val = " << unsigned(sum) << "\n";
    for (int_fast16_t i = 0; i < block_size; ++i) {
        out_block[i] = dc_val;
    }
}

inline void kvz_intra_pred_planar_generic(
    const int_fast8_t log2_width,
    const kvz_pixel* const ref_top,
    const kvz_pixel* const ref_left,
    kvz_pixel* const dst)
{
    assert(log2_width >= 2 && log2_width <= 5);

    const int_fast8_t width = 1 << log2_width;
    const kvz_pixel top_right = ref_top[width + 1];
    const kvz_pixel bottom_left = ref_left[width + 1];

#if 1
    // Unoptimized version for reference.
    for (int y = 0; y < width; ++y) {
        for (int x = 0; x < width; ++x) {
            int_fast16_t hor = (width - 1 - x) * ref_left[y + 1] + (x + 1) * top_right;
            int_fast16_t ver = (width - 1 - y) * ref_top[x + 1] + (y + 1) * bottom_left;
            dst[y * width + x] = (ver + hor + width) >> (log2_width + 1);
        }
    }
#else
    int_fast16_t top[32];
    for (int i = 0; i < width; ++i) {
        top[i] = ref_top[i + 1] << log2_width;
    }

    for (int y = 0; y < width; ++y) {
        int_fast16_t hor = (ref_left[y + 1] << log2_width) + width;
        for (int x = 0; x < width; ++x) {
            hor += top_right - ref_left[y + 1];
            top[x] += bottom_left - ref_top[x + 1];
            dst[y * width + x] = (hor + top[x]) >> (log2_width + 1);
        }
    }
#endif
}

inline void kvz_intra_pred_filtered_dc_generic(
    const int_fast8_t log2_width,
    const kvz_pixel* const ref_top,
    const kvz_pixel* const ref_left,
    kvz_pixel* const out_block)
{
    assert(log2_width >= 2 && log2_width <= 5);

    const int_fast8_t width = 1 << log2_width;

    int_fast16_t sum = 0;
    for (int_fast8_t i = 0; i < width; ++i) {
        sum += ref_top[i + 1];
        sum += ref_left[i + 1];
    }

    const kvz_pixel dc_val = (sum + width) >> (log2_width + 1);

    // Filter top-left with ([1 2 1] / 4)
    out_block[0] = (ref_left[1] + 2 * dc_val + ref_top[1] + 2) / 4;

    // Filter rest of the boundary with ([1 3] / 4)
    for (int_fast8_t x = 1; x < width; ++x) {
        out_block[x] = (ref_top[x + 1] + 3 * dc_val + 2) / 4;
    }
    for (int_fast8_t y = 1; y < width; ++y) {
        out_block[y * width] = (ref_left[y + 1] + 3 * dc_val + 2) / 4;
        for (int_fast8_t x = 1; x < width; ++x) {
            out_block[y * width + x] = dc_val;
        }
    }
}



inline void kvz_angular_pred_generic(
    const int_fast8_t log2_width,
    const int_fast8_t intra_mode,
    const kvz_pixel* const in_ref_above,
    const kvz_pixel* const in_ref_left,
    kvz_pixel* const dst)
{
    assert(log2_width >= 2 && log2_width <= 5);
    assert(intra_mode >= 2 && intra_mode <= 34);

    static const int8_t modedisp2sampledisp[9] = { 0, 2, 5, 9, 13, 17, 21, 26, 32 };
    static const int16_t modedisp2invsampledisp[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / sampledisp

                                                      // Temporary buffer for modes 11-25.
                                                      // It only needs to be big enough to hold indices from -width to width-1.
    kvz_pixel tmp_ref[2 * 32];
    const int_fast8_t width = 1 << log2_width;

    // Whether to swap references to always project on the left reference row.
    const bool vertical_mode = intra_mode >= 18;
    // Modes distance to horizontal or vertical mode.
    const int_fast8_t mode_disp = vertical_mode ? intra_mode - 26 : 10 - intra_mode;
    // Sample displacement per column in fractions of 32.
    const int_fast8_t sample_disp = (mode_disp < 0 ? -1 : 1) * modedisp2sampledisp[abs(mode_disp)];

    // Pointer for the reference we are interpolating from.
    const kvz_pixel* ref_main;
    // Pointer for the other reference.
    const kvz_pixel* ref_side;

    // Set ref_main and ref_side such that, when indexed with 0, they point to
    // index 0 in block coordinates.
    if (sample_disp < 0) {
        // Negative sample_disp means, we need to use both references.

        ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
        ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;

        // Move the reference pixels to start from the middle to the later half of
        // the tmp_ref, so there is room for negative indices.
        for (int_fast8_t x = -1; x < width; ++x) {
            tmp_ref[x + width] = ref_main[x];
        }
        // Get a pointer to block index 0 in tmp_ref.
        ref_main = &tmp_ref[width];

        // Extend the side reference to the negative indices of main reference.
        int_fast32_t col_sample_disp = 128; // rounding for the ">> 8"
        int_fast16_t inv_abs_sample_disp = modedisp2invsampledisp[abs(mode_disp)];
        int_fast8_t most_negative_index = (width * sample_disp) >> 5;
        for (int_fast8_t x = -2; x >= most_negative_index; --x) {
            col_sample_disp += inv_abs_sample_disp;
            int_fast8_t side_index = col_sample_disp >> 8;
            tmp_ref[x + width] = ref_side[side_index - 1];
        }
    }
    else {
        // sample_disp >= 0 means we don't need to refer to negative indices,
        // which means we can just use the references as is.
        ref_main = (vertical_mode ? in_ref_above : in_ref_left) + 1;
        ref_side = (vertical_mode ? in_ref_left : in_ref_above) + 1;
    }

    if (sample_disp != 0) {
        // The mode is not horizontal or vertical, we have to do interpolation.

        int_fast16_t delta_pos = 0;
        for (int_fast8_t y = 0; y < width; ++y) {
            delta_pos += sample_disp;
            int_fast8_t delta_int = delta_pos >> 5;
            int_fast8_t delta_fract = delta_pos & (32 - 1);

            if (delta_fract) {
                // Do linear filtering
                for (int_fast8_t x = 0; x < width; ++x) {
                    kvz_pixel ref1 = ref_main[x + delta_int];
                    kvz_pixel ref2 = ref_main[x + delta_int + 1];
                    dst[y * width + x] = ((32 - delta_fract) * ref1 + delta_fract * ref2 + 16) >> 5;
                }
            }
            else {
                // Just copy the integer samples
                for (int_fast8_t x = 0; x < width; x++) {
                    dst[y * width + x] = ref_main[x + delta_int];
                }
            }
        }
    }
    else {
        // Mode is horizontal or vertical, just copy the pixels.

        for (int_fast8_t y = 0; y < width; ++y) {
            for (int_fast8_t x = 0; x < width; ++x) {
                dst[y * width + x] = ref_main[x];
            }
        }
    }

    // Flip the block if this is was a horizontal mode.
    if (!vertical_mode) {
        for (int_fast8_t y = 0; y < width - 1; ++y) {
            for (int_fast8_t x = y + 1; x < width; ++x) {
                SWAP(dst[y * width + x], dst[x * width + y], kvz_pixel);
            }
        }
    }
}


inline void intra_post_process_angular(
    unsigned width,
    unsigned stride,
    const kvz_pixel* ref,
    kvz_pixel* block)
{
    kvz_pixel ref2 = ref[0];
    for (unsigned i = 0; i < width; i++) {
        kvz_pixel val = block[i * stride];
        kvz_pixel ref1 = ref[i + 1];
        block[i * stride] = CLIP_TO_PIXEL(val + ((ref1 - ref2) >> 1));
    }
}


inline void kvz_intra_predict(
    kvz_intra_references* refs,
    int_fast8_t log2_width,
    int_fast8_t mode,
    color_t color,
    kvz_pixel* dst,
    bool filter_boundary)
{
    const int_fast8_t width = 1 << log2_width;
    /*-----------------------------------------------------*/
    const kvz_intra_ref* used_ref = &refs->ref;
    if (color != COLOR_Y || mode == 1 || width == 4) {
        // For chroma, DC and 4x4 blocks, always use unfiltered reference.
    }
    else if (mode == 0) {
        // Otherwise, use filtered for planar.
        used_ref = &refs->filtered_ref;
    }
    else {
        // Angular modes use smoothed reference pixels, unless the mode is close
        // to being either vertical or horizontal.
        static const int kvz_intra_hor_ver_dist_thres[5] = { 0, 7, 1, 0, 0 };
        int filter_threshold = kvz_intra_hor_ver_dist_thres[kvz_math_floor_log2(width) - 2];
        int dist_from_vert_or_hor = MIN(abs(mode - 26), abs(mode - 10));
        if (dist_from_vert_or_hor > filter_threshold) {
            used_ref = &refs->filtered_ref;
        }
    }
    /*--------------------------------------------------------*/
    if (used_ref == &refs->filtered_ref && !refs->filtered_initialized) {
        intra_filter_reference(log2_width, refs);
    }
    /*--------------------------------------------------------*/
    if (mode == 0) {
        kvz_intra_pred_planar_generic(log2_width, used_ref->top, used_ref->left, dst);
    }
    else if (mode == 1) {
        // Do extra post filtering for edge pixels of luma DC mode.
        if (color == COLOR_Y && width < 32) {
            kvz_intra_pred_filtered_dc_generic(log2_width, used_ref->top, used_ref->left, dst);
        }
        else {
            intra_pred_dc(log2_width, used_ref->top, used_ref->left, dst);
        }
    }
    else {
        kvz_angular_pred_generic(log2_width, mode, used_ref->top, used_ref->left, dst);
        if (color == COLOR_Y && width < 32 && filter_boundary) {
            if (mode == 10) {
                intra_post_process_angular(width, 1, used_ref->top, dst);
            }
            else if (mode == 26) {
                intra_post_process_angular(width, width, used_ref->left, dst);
            }
        }
    }
}