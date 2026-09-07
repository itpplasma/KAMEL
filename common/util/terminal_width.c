#if defined(__linux__) || defined(__APPLE__)
#include <sys/ioctl.h>
#include <unistd.h>
#endif

int kamel_get_terminal_width_for_fd(int const file_descriptor) {
#if defined(__linux__) || defined(__APPLE__)
    struct winsize window = {0};

    int const ret = ioctl(file_descriptor, TIOCGWINSZ, &window);

    if (file_descriptor < 0 || ret == -1 || window.ws_col == 0) {
        return 0;
    }

    return (int)window.ws_col;
#else
    (void)file_descriptor;
    return 0;
#endif
}

int kamel_get_terminal_width(void) {
#if defined(__linux__) || defined(__APPLE__)
    return kamel_get_terminal_width_for_fd(STDOUT_FILENO);
#else
    return 0;
#endif
}
