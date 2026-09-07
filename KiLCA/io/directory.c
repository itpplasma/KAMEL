/* Keep the platform-specific struct dirent layout behind its native C API. */
#include <dirent.h>
#include <errno.h>
#include <stddef.h>
#include <string.h>

int kilca_readdir_name(void* directory, char* name, size_t capacity) {
    errno = 0;
    struct dirent* entry = readdir(directory);
    if (entry == NULL)
        return errno == 0 ? 0 : -1;
    size_t length = strlen(entry->d_name);
    if (length >= capacity)
        return -1;
    memcpy(name, entry->d_name, length + 1);
    return 1;
}
