#ifndef BGEN_PROGRESS_H
#define BGEN_PROGRESS_H

/* Minimal, dependency-free progress bar: percentage, rate-derived ETA, no
 * external library. Replaces almosthere (athr), whose CMake packaging kept
 * breaking downstream consumers with dangling internal link targets. */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#if defined(_WIN32)
#include <io.h>
#define bgen_progress_isatty(fd) _isatty(fd)
#define bgen_progress_fileno(f) _fileno(f)
#else
#include <unistd.h>
#define bgen_progress_isatty(fd) isatty(fd)
#define bgen_progress_fileno(f) fileno(f)
#endif

struct bgen_progress
{
    char const* label;
    uint64_t    total;
    uint64_t    done;
    time_t      start;
    time_t      last_render;
    int         is_tty;
};

static void bgen_progress_render(struct bgen_progress const* p, time_t now)
{
    double elapsed = (double)(now - p->start);
    double frac    = p->total ? (double)p->done / (double)p->total : 1.0;
    if (frac > 1.0)
        frac = 1.0;

    char eta[16] = "--:--:--";
    if (p->done > 0 && p->done < p->total) {
        double rate = (double)p->done / (elapsed > 0.0 ? elapsed : 1.0);
        if (rate > 0.0) {
            long secs = (long)((double)(p->total - p->done) / rate);
            snprintf(eta, sizeof(eta), "%02ld:%02ld:%02ld", secs / 3600, (secs / 60) % 60,
                     secs % 60);
        }
    } else if (p->done >= p->total) {
        snprintf(eta, sizeof(eta), "00:00:00");
    }

    int const width  = 24;
    int const filled = (int)(frac * width);
    char      bar[width + 1];
    int       i = 0;
    for (; i < filled; ++i)
        bar[i] = '#';
    for (; i < width; ++i)
        bar[i] = '-';
    bar[width] = '\0';

    fprintf(stderr, "%s%s [%s] %5.1f%% (%" PRIu64 "/%" PRIu64 ") ETA %s%s",
            p->is_tty ? "\r" : "", p->label, bar, frac * 100.0, p->done, p->total, eta,
            p->is_tty ? "" : "\n");
    fflush(stderr);
}

static void bgen_progress_start(struct bgen_progress* p, uint64_t total, char const* label)
{
    p->label       = label;
    p->total       = total;
    p->done        = 0;
    p->start       = time(NULL);
    p->last_render = p->start;
    p->is_tty      = bgen_progress_isatty(bgen_progress_fileno(stderr));
    bgen_progress_render(p, p->start);
}

static void bgen_progress_update(struct bgen_progress* p, uint64_t n)
{
    p->done += n;
    time_t now = time(NULL);
    /* Throttled to once a second, so a fast scan doesn't spend more time
     * rendering than working. */
    if (now != p->last_render || p->done >= p->total) {
        p->last_render = now;
        bgen_progress_render(p, now);
    }
}

static void bgen_progress_finish(struct bgen_progress* p)
{
    p->done = p->total;
    bgen_progress_render(p, time(NULL));
    fputc('\n', stderr);
}

#endif
