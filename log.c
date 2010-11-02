#include <stdarg.h>
#include "sting.h"
#include "log.h"

int open_log(const char *fname)
{
    logfp = fopen(fname, "a+");

    if (!logfp)
    {
        fprintf(stderr, "Error: Unable to open log file\n");
        return 1;
    }

    return 0;
}

int flog(const char *fmt, ...)
{
    va_list ap;

    if (logfp != NULL)
    {
        va_start(ap, fmt);
        vfprintf(logfp, fmt, ap);
        va_end(ap);
    }

    if (echolog)
    {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        va_end(ap);
    }

    return 0;
}
