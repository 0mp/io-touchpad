/*
 * Based on the evtest source code (version 1.32 downloaded from Debian).
 */

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#define _GNU_SOURCE /* for asprintf */
#include <stdio.h>
#include <stdint.h>

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <linux/version.h>
#include <linux/input.h>

#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <getopt.h>
#include <ctype.h>

#include "touchpadlib.h"

#define BITS_PER_LONG (sizeof(long) * 8)
#define NBITS(x) ((((x)-1)/BITS_PER_LONG)+1)
#define OFF(x)  ((x)%BITS_PER_LONG)
#define BIT(x)  (1UL<<OFF(x))
#define LONG(x) ((x)/BITS_PER_LONG)
#define test_bit(bit, array)    ((array[LONG(bit)] >> OFF(bit)) & 1)

#define DEV_INPUT_EVENT "/dev/input"
#define EVENT_DEV_NAME "event"

#ifndef EV_SYN
#define EV_SYN 0
#endif

#define NAME_ELEMENT(element) [element] = #element

static const char * const events[EV_MAX + 1] = {
    [0 ... EV_MAX] = NULL,
    NAME_ELEMENT(EV_ABS),
};

static const int maxval[EV_MAX + 1] = {
    [0 ... EV_MAX] = -1,
    [EV_ABS] = ABS_MAX,
};

static const char * const absval[6] = {
    "Value", "Min  ",
    "Max  ", "Fuzz ",
    "Flat ", "Resolution "
};

static const char * const absolutes[ABS_MAX + 1] = {
    [0 ... ABS_MAX] = NULL,
    NAME_ELEMENT(ABS_X),            NAME_ELEMENT(ABS_Y),
    NAME_ELEMENT(ABS_Z),            NAME_ELEMENT(ABS_RX),
    NAME_ELEMENT(ABS_RY),           NAME_ELEMENT(ABS_RZ),
    NAME_ELEMENT(ABS_THROTTLE),     NAME_ELEMENT(ABS_RUDDER),
    NAME_ELEMENT(ABS_WHEEL),        NAME_ELEMENT(ABS_GAS),
    NAME_ELEMENT(ABS_BRAKE),        NAME_ELEMENT(ABS_HAT0X),
    NAME_ELEMENT(ABS_HAT0Y),        NAME_ELEMENT(ABS_HAT1X),
    NAME_ELEMENT(ABS_HAT1Y),        NAME_ELEMENT(ABS_HAT2X),
    NAME_ELEMENT(ABS_HAT2Y),        NAME_ELEMENT(ABS_HAT3X),
    NAME_ELEMENT(ABS_HAT3Y),        NAME_ELEMENT(ABS_PRESSURE),
    NAME_ELEMENT(ABS_DISTANCE),     NAME_ELEMENT(ABS_TILT_X),
    NAME_ELEMENT(ABS_TILT_Y),       NAME_ELEMENT(ABS_TOOL_WIDTH),
    NAME_ELEMENT(ABS_VOLUME),       NAME_ELEMENT(ABS_MISC),
#ifdef ABS_MT_BLOB_ID
    NAME_ELEMENT(ABS_MT_TOUCH_MAJOR),
    NAME_ELEMENT(ABS_MT_TOUCH_MINOR),
    NAME_ELEMENT(ABS_MT_WIDTH_MAJOR),
    NAME_ELEMENT(ABS_MT_WIDTH_MINOR),
    NAME_ELEMENT(ABS_MT_ORIENTATION),
    NAME_ELEMENT(ABS_MT_POSITION_X),
    NAME_ELEMENT(ABS_MT_POSITION_Y),
    NAME_ELEMENT(ABS_MT_TOOL_TYPE),
    NAME_ELEMENT(ABS_MT_BLOB_ID),
#endif
#ifdef ABS_MT_TRACKING_ID
    NAME_ELEMENT(ABS_MT_TRACKING_ID),
#endif
#ifdef ABS_MT_PRESSURE
    NAME_ELEMENT(ABS_MT_PRESSURE),
#endif
#ifdef ABS_MT_SLOT
    NAME_ELEMENT(ABS_MT_SLOT),
#endif
#ifdef ABS_MT_TOOL_X
    NAME_ELEMENT(ABS_MT_TOOL_X),
    NAME_ELEMENT(ABS_MT_TOOL_Y),
    NAME_ELEMENT(ABS_MT_DISTANCE),
#endif
};

static const char * const * const names[EV_MAX + 1] = {
    [0 ... EV_MAX] = NULL,
    [EV_SYN] = events,
    [EV_ABS] = absolutes,
};

/**
 * Filter for the AutoDevProbe scandir on /dev/input.
 *
 * @param dir The current directory entry provided by scandir.
 *
 * @return Non-zero if the given directory entry starts with "event", or zero
 * otherwise.
 */
static int is_event_device(const struct dirent *dir) {
    return strncmp(EVENT_DEV_NAME, dir->d_name, 5) == 0;
}

/**
 * Scans all /dev/input/event*, display them and ask the user which one to
 * open.
 *
 * @return The event device file name of the device file selected. This
 * string is allocated and must be freed by the caller.
 */
static char* scan_devices(void)
{
    struct dirent **namelist;
    int i, ndev, devnum;
    char *filename;

    ndev = scandir(DEV_INPUT_EVENT, &namelist, is_event_device, versionsort);
    if (ndev <= 0)
        return NULL;

    fprintf(stderr, "Available devices:\n");

    for (i = 0; i < ndev; i++)
    {
        char fname[64];
        int fd = -1;
        char name[256] = "???";

        snprintf(fname, sizeof(fname),
             "%s/%s", DEV_INPUT_EVENT, namelist[i]->d_name);
        fd = open(fname, O_RDONLY);
        if (fd < 0)
            continue;
        ioctl(fd, EVIOCGNAME(sizeof(name)), name);

        // @TODO We can try to detect touchpad device here if we examine name for the 'touchpad' substring.

        fprintf(stderr, "%s:    %s\n", fname, name);
        close(fd);
        free(namelist[i]);
    }

    fprintf(stderr, "Select the device event number [0-%d]: ", ndev - 1);
    scanf("%d", &devnum);

    if (devnum >= ndev || devnum < 0)
        return NULL;

    asprintf(&filename, "%s/%s%d",
         DEV_INPUT_EVENT, EVENT_DEV_NAME,
         devnum);

    return filename;
}

static void print_absdata(int fd, int axis)
{
    int abs[6] = {0};
    int k;

    ioctl(fd, EVIOCGABS(axis), abs);
    for (k = 0; k < 6; k++)
        if ((k < 3) || abs[k])
            printf("      %s %6d\n", absval[k], abs[k]);
}

static void print_repdata(int fd)
{
    int i;
    unsigned int rep[2];

    ioctl(fd, EVIOCGREP, rep);

    for (i = 0; i <= REP_MAX; i++) {
        printf("    Repeat code %d (%s)\n", i, names[EV_REP] ? (names[EV_REP][i] ? names[EV_REP][i] : "?") : "?");
        printf("      Value %6d\n", rep[i]);
    }

}

static inline const char* typename(unsigned int type)
{
    return (type <= EV_MAX && events[type]) ? events[type] : "?";
}

static inline const char* codename(unsigned int type, unsigned int code)
{
    return (type <= EV_MAX && code <= maxval[type] && names[type] && names[type][code]) ? names[type][code] : "?";
}

/**
 * Print static device information (no events). This information includes
 * version numbers, device name and all bits supported by this device.
 *
 * @param fd The file descriptor to the device.
 * @return 0 on success or 1 otherwise.
 */
static int print_device_info(int fd)
{
    unsigned int type, code;
    int version;
    unsigned short id[4];
    char name[256] = "Unknown";
    unsigned long bit[EV_MAX][NBITS(KEY_MAX)];

    if (ioctl(fd, EVIOCGVERSION, &version)) {
        perror("ERROR: touchpadlib: can't get version");
        return 1;
    }

    printf("Input driver version is %d.%d.%d\n",
        version >> 16, (version >> 8) & 0xff, version & 0xff);

    ioctl(fd, EVIOCGID, id);
    printf("Input device ID: bus 0x%x vendor 0x%x product 0x%x version 0x%x\n",
        id[ID_BUS], id[ID_VENDOR], id[ID_PRODUCT], id[ID_VERSION]);

    ioctl(fd, EVIOCGNAME(sizeof(name)), name);
    printf("Input device name: \"%s\"\n", name);

    memset(bit, 0, sizeof(bit));
    ioctl(fd, EVIOCGBIT(0, EV_MAX), bit[0]);
    printf("Supported events:\n");

    for (type = 0; type < EV_MAX; type++) {
        if (test_bit(type, bit[0]) && type != EV_REP) {
            printf("  Event type %d (%s)\n", type, typename(type));
            if (type == EV_SYN) continue;
            ioctl(fd, EVIOCGBIT(type, KEY_MAX), bit[type]);
            for (code = 0; code < KEY_MAX; code++)
                if (test_bit(code, bit[type])) {
                    printf("    Event code %d (%s)\n", code, codename(type, code));
                    if (type == EV_ABS)
                        print_absdata(fd, code);
                }
        }
    }

    if (test_bit(EV_REP, bit[0])) {
        printf("Key repeat handling:\n");
        printf("  Repeat type %d (%s)\n", EV_REP, events[EV_REP] ?  events[EV_REP] : "?");
        print_repdata(fd);
    }
    return 0;
}

/**
 * Reset fields of a touchpad event.
 *
 * @param event Pointer to the touchpad_event you want to reset.
 */
static void reset_touchpad_event(struct touchpad_event *event)
{
    event->x = -1;
    event->y = -1;
    event->pressure = -1;
    event->seconds = 0;
    event->useconds = 0;
}

/**
 * Print the content of struct touchpad_event in a human readable format.
 *
 * @param event A pointer to the touchpad_event which we want to print.
 */
void print_event(const struct touchpad_event *event)
{
    printf("ABS_X %d\tABS_Y %d\tABS_PRESSURE %d\t"
            "seconds %ld\tmiliseconds %ld\n",
            event->x, event->y, event->pressure,
            event->seconds, event->useconds);
}

/**
 * Get the next device input event as it comes in.
 *
 * @param fd The file descriptor to the device.
 * @param touchpad_event The structure where the data will be saved.
 * @return 0 on success or 1 otherwise.
 */
int fetch_touchpad_event(int fd, struct touchpad_event *event)
{
    struct input_event ev[64];
    int i, rd;

    rd = read(fd, ev, sizeof(struct input_event) * 64);

    if (rd < (int) sizeof(struct input_event)) {
        fprintf(stderr, "ERROR: Expected %d bytes, got %d.\n",
                (int) sizeof(struct input_event), rd);
        perror("\nERROR: touchpadlib: error reading.");
        return 1;
    }

    reset_touchpad_event(event);

    for (i = 0; i < rd / sizeof(struct input_event); i++) {
        unsigned int type, code;

        type = ev[i].type;
        code = ev[i].code;

        // Set seconds only during the first iteration over read events.
        if (i == 0) {
            event->seconds = ev[i].time.tv_sec;
            event->useconds = ev[i].time.tv_usec;
        }

        if (type == EV_ABS) {
            if (code == ABS_X) {
                event->x = ev[i].value;
            } else if (code == ABS_Y) {
                event->y = ev[i].value;
            } else if (code == ABS_PRESSURE) {
                event->pressure = ev[i].value;
            }
        }
    }
    return 0;
}

/**
 * Grab and immediately ungrab the device.
 *
 * @param fd The file descriptor to the device.
 * @return 0 if the grab was successful, or 1 otherwise.
 */
static int test_grab(int fd)
{
    int rc;

    rc = ioctl(fd, EVIOCGRAB, (void*)1);

    if (!rc)
        ioctl(fd, EVIOCGRAB, (void*)0);

    return rc;
}

/**
 * @TODO
 * Enter capture mode. The requested event device will be monitored, and any
 * captured events will be decoded and printed on the console.
 * Initialized usage of the touchpadlib library.
 *
 * The function will:
 * - Check if the user has privileges of root.
 * - Show the list of available device so that the user can choose the touchpad
 *   device. (@TODO It could be done automaically. See issue #4)
 * - Open the device for reading.
 *
 * @return File descriptor to the device selected by the user.
 *  -1 if the function fails to return a proper file descriptor.
 */
int initalize_touchpadlib_usage()
{
    int fd;
    char *filename = NULL;

    if (!has_root_privileges())
        fprintf(stderr, "WARNING: Not running as root, no devices may be available.\n");

    filename = scan_devices();

    if (!filename)
        return EXIT_FAILURE;

    if ((fd = open(filename, O_RDONLY)) < 0) {
        perror("ERROR: touchpadlib");
        if (errno == EACCES && getuid() != 0)
            fprintf(stderr, "ERROR: You do not have access to %s. Try "
                    "running as root instead.\n",
                    filename);
        goto error;
    }

    if (!isatty(fileno(stdout)))
        setbuf(stdout, NULL);

    // if (print_device_info(fd))
    //     goto error;

    if (test_grab(fd))
        goto error;

    free(filename);

    return fd;

error:
    free(filename);
    return -1;
}

/**
 * Check if the program is ran by root.
 *
 * @return 1 if run by run, 0 otherwise.
 */
int has_root_privileges(void)
{
    return getuid() == 0;
}

/**
 * Allocate a new struct touchpad_event.
 *
 * It has to be freed by the user later on.
 *
 * @return Pointer to the allocated struct. NULL if malloc failed.
 */
struct touchpad_event *new_event() {
    return malloc(sizeof(struct touchpad_event));
}

/**
 * Free the touchpad_event.
 *
 * @param event The event to be freed.
 */
void erase_event(struct touchpad_event *event) {
    free(event);
}

/**
 * Getter for the coordinate x within a touchpad_event.
 *
 * @param event Pointer to the event.
 * @return The x coordinate.
 */
int get_x(struct touchpad_event *event) {
    return event->x;
}

/**
 * Getter for the coordinate y within a touchpad_event.
 *
 * @param event Pointer to the event.
 * @return The y coordinate.
 */
int get_y(struct touchpad_event *event) {
    return event->y;
}

/**
 * Getter for the pressure within a touchpad_event.
 *
 * @param event Pointer to the event.
 * @return The value of the pressure field.
 */
int get_pressure(struct touchpad_event *event) {
    return event->pressure;
}

/**
 * Getter for the seconds within a touchpad_event.
 *
 * @param event Pointer to the event.
 * @return The value of the seconds field.
 */
int get_seconds(struct touchpad_event *event) {
    return event->seconds;
}

/**
 * Getter for the miliseconds within a touchpad_event.
 *
 * @param event Pointer to the event.
 * @return The value of the miliseconds field.
 */
int get_useconds(struct touchpad_event *event) {
    return event->useconds;
}
