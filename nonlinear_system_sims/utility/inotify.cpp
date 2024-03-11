// Auto detect file modification
// for shader auto-reloading
#include "inotify.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <limits.h>
#include <unistd.h>

enum{
    MAX_EVENTS = 64, /*Max. number of events to process at one go*/
    LEN_NAME   = 16, /*Assuming that the length of the filename won't exceed 16 bytes*/
    EVENT_SIZE = ( sizeof (struct inotify_event) ), /*size of one event*/
    BUF_LEN    = ( MAX_EVENTS * ( EVENT_SIZE + LEN_NAME )) /*buffer to store the data of events*/
};

Inotify* Inotify::m_singleton = 0;

int isready(int fd, int usec)
{
int rc;
fd_set fds;
struct timeval tv;
FD_ZERO(&fds);
FD_SET(fd,&fds);
tv.tv_sec = 0;tv.tv_usec = usec;
rc = select(fd+1, &fds, NULL, NULL, &tv);
if (rc < 0) //error
return -1;
return FD_ISSET(fd,&fds) ? 1 : 0;
}

Inotify::~Inotify(){
    // Clean up
    if(_buffer) delete [] _buffer;
    inotify_rm_watch( _fd, _wd );
    close( _fd );
}

Inotify::Inotify(const char* path):_buffer(nullptr),_fd(-1),_wd(-1){
    // Initialize Inotify
#ifdef IN_NONBLOCK
    _fd = inotify_init1( IN_NONBLOCK );
#else
    _fd = inotify_init();
#endif

    if ( _fd < 0 ) {
        printf("Couldn't initialize non-blocking inotify\n");
    }

    // add watch to starting directory
    _wd = inotify_add_watch(_fd, path, IN_CREATE | IN_MODIFY | IN_DELETE);

    if (_wd == -1)
    {
        printf("Couldn't add watch to %s\n",path);
    }
    else
    {
        printf("Watching:: %s\n",path);
    }
    if(!_buffer) _buffer = new char[BUF_LEN];
}

bool Inotify::_update(){

    int length = 0, i = 0;

    if(isready(_fd, 1)) length = read( _fd, _buffer, BUF_LEN );

    if ( length < 0 ) {
        perror("read");
    }

    while ( i < length ) {
        //perror("is reading events...");
        struct inotify_event *event = ( struct inotify_event * ) &_buffer[ i ];
        if ( event->len ) {
            if ( event->mask & IN_CREATE) {
                if (event->mask & IN_ISDIR)
                {/*printf( "The directory %s was Created.\n", event->name );*/return true;}
                else
                {/*printf( "The file %s was Created with WD %d\n", event->name, event->wd );*/return true;}
            }

            if ( event->mask & IN_MODIFY) {
                if (event->mask & IN_ISDIR)
                {/*printf( "The directory %s was modified.\n", event->name );*/return true;}
                else
                {/*printf( "The file %s was modified with WD %d\n", event->name, event->wd );*/return true;}
            }

            if ( event->mask & IN_DELETE) {
                if (event->mask & IN_ISDIR)
                {/*printf( "The directory %s was deleted.\n", event->name );*/return true;}
                else
                {/*printf( "The file %s was deleted with WD %d\n", event->name, event->wd );*/return true;}
            }

            i += EVENT_SIZE + event->len;
        }
    }

    return false;
}
