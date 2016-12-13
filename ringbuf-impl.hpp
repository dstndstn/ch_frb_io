#include <iostream>
#include <deque>

template <class T>
Ringbuf<T>::Ringbuf(int maxsize) :
    _deleter(this),
    _q(),
    _live(0),
    _maxsize(maxsize)
{
    pthread_mutex_init(&this->_q_lock, NULL);
    pthread_mutex_init(&this->_live_lock, NULL);
}

template <class T>
Ringbuf<T>::~Ringbuf() {
    pthread_mutex_destroy(&this->_q_lock);
    pthread_mutex_destroy(&this->_live_lock);
}

template <class T>
std::shared_ptr<T> Ringbuf<T>::push(T* t) {
    // Is there room?  (This tries to make room if not...)
    bool can = can_push();
    if (!can)
        return std::shared_ptr<T>();
    //std::cout << "Creating shared_ptr..." << std::endl;
    std::shared_ptr<T> p(t, _deleter);
    pthread_mutex_lock(&this->_live_lock);
    _live++;
    pthread_mutex_unlock(&this->_live_lock);
    pthread_mutex_lock(&this->_q_lock);
    _q.push_back(p);
    pthread_mutex_unlock(&this->_q_lock);
    //std::cout << "Now " << _live << " objects are live" << std::endl;
    return p;
}

template <class T>
std::shared_ptr<T> Ringbuf<T>::pop() {
    //std::cout << "Popping..." << std::endl;
    std::shared_ptr<T> p;
    pthread_mutex_lock(&this->_q_lock);
    if (_q.empty()) {
        //std::cout << "Pop: empty" << std::endl;
    } else {
        p = _q.pop_front();
        //std::cout << "Pop: returning " << *p << std::endl;
    }
    pthread_mutex_unlock(&this->_q_lock);
    return p;
}

template <class T>
std::vector<std::shared_ptr<T> > Ringbuf<T>::snapshot(bool (*testFunc)(const std::shared_ptr<T>)) {
    std::vector<std::shared_ptr<T> > vec;
    snapshot(vec, testFunc);
    return vec;
}

template <class T>
void Ringbuf<T>::snapshot(std::vector<std::shared_ptr<T> > &vec,
                          bool (*testFunc)(const std::shared_ptr<T>)) {
    pthread_mutex_lock(&this->_q_lock);
    for (auto it = _q.begin(); it != _q.end(); it++) {
        if (!testFunc || testFunc(*it)) {
	    vec.push_back(*it);
        }
    }
    pthread_mutex_unlock(&this->_q_lock);
}

template <class T>
void Ringbuf<T>::dropping(std::shared_ptr<T> t) {}

// Called by the RingbufDeleter when a shared_ptr is deleted
template <class T>
void Ringbuf<T>::deleted(T* t) {
    //std::cout << "Deleting object: " << *t << std::endl;
    pthread_mutex_lock(&this->_live_lock);
    _live--;
    pthread_mutex_unlock(&this->_live_lock);
    //std::cout << "Now " << _live << " objects are live" << std::endl;
    delete t;
}

template <class T>
bool Ringbuf<T>::can_push() {
    while (1) {
        pthread_mutex_lock(&this->_live_lock);
        size_t nlive = _live;
        pthread_mutex_unlock(&this->_live_lock);
        if (nlive < _maxsize)
            break;
        //std::cout << "Push: _live >= _maxsize." << " (" << nlive << " >= " << _maxsize << ")" << std::endl;
        pthread_mutex_lock(&this->_q_lock);

        if (_q.empty()) {
            //std::cout << "Ring buffer empty but still too many live elements -- push fails" << std::endl;
            pthread_mutex_unlock(&this->_q_lock);
            return false;
        }
        //std::cout << "Dropping an element..." << std::endl;
        std::shared_ptr<T> p = _q.front();
        _q.pop_front();
        pthread_mutex_unlock(&this->_q_lock);
        dropping(p);
        p.reset();
        //std::cout << "Now " << _live << " live" << std::endl;
    }
    return true;
}

