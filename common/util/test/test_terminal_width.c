#define _XOPEN_SOURCE 600

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

int kamel_get_terminal_width_for_fd(int file_descriptor);

static int check_width(int file_descriptor, int expected_width)
{
    int actual_width = kamel_get_terminal_width_for_fd(file_descriptor);

    if (actual_width != expected_width) {
        fprintf(stderr, "expected terminal width %d, got %d\n", expected_width, actual_width);
        return 0;
    }

    return 1;
}

int main(void)
{
    int master_fd;
    int slave_fd;
    char *slave_name;
    struct winsize window = {0};

    if (!check_width(-1, 0)) {
        return EXIT_FAILURE;
    }

    master_fd = posix_openpt(O_RDWR | O_NOCTTY);
    if (master_fd == -1) {
        return EXIT_SUCCESS;
    }

    if (grantpt(master_fd) == -1 || unlockpt(master_fd) == -1) {
        close(master_fd);
        return EXIT_SUCCESS;
    }

    slave_name = ptsname(master_fd);
    if (slave_name == NULL) {
        close(master_fd);
        return EXIT_SUCCESS;
    }

    slave_fd = open(slave_name, O_RDWR | O_NOCTTY);
    if (slave_fd == -1) {
        close(master_fd);
        return EXIT_SUCCESS;
    }

    window.ws_col = 80;
    if (ioctl(slave_fd, TIOCSWINSZ, &window) == -1 || !check_width(slave_fd, 80)) {
        close(slave_fd);
        close(master_fd);
        return EXIT_FAILURE;
    }

    window.ws_col = 120;
    if (ioctl(slave_fd, TIOCSWINSZ, &window) == -1 || !check_width(slave_fd, 120)) {
        close(slave_fd);
        close(master_fd);
        return EXIT_FAILURE;
    }

    close(slave_fd);
    close(master_fd);
    return EXIT_SUCCESS;
}
