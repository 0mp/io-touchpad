#ifndef TOUCHPADLIB_H
#define TOUCHPADLIB_H

#include <sys/time.h>
#include <linux/input.h>

struct touchpad_event {
    int32_t     x;          // The absolute coordinate x.
    int32_t     y;          // The absolute coordinate y.
    int32_t     pressure;   // Pressure applied to touchpad surface.
    time_t      seconds;    // Timestamp of the event occurance in seconds.
    suseconds_t useconds;   // Additional precision for in miliseconds.
};

/* struct touchpad_event interface. */
void erase_event(struct touchpad_event *event);
int get_x(struct touchpad_event *event);
int get_y(struct touchpad_event *event);
int get_pressure(struct touchpad_event *event);
int get_seconds(struct touchpad_event *event);
int get_useconds(struct touchpad_event *event);
struct touchpad_event *new_event(void);

/* touchpadlib interface */
int fetch_touchpad_event(int fd, struct touchpad_event *event);
int has_root_privileges(void);
int initalize_touchpadlib_usage(void);
void print_event(const struct touchpad_event *event);

#endif /* TOUCHPADLIB_H */
