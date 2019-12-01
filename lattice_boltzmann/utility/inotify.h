#ifndef INOTIFY_H
#define INOTIFY_H

class Inotify
{
public:
    static void Init(const char* path){m_singleton = new Inotify(path);}
    static bool Update() {return m_singleton->_update();}
    static void Finalize(){delete m_singleton;}
private:
    static Inotify* m_singleton;
    int _fd, _wd;
    char* _buffer;

    Inotify(const char*);
    ~Inotify();
    bool _update();
};

#endif
