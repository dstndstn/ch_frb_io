#include <iostream>
#include <deque>

#include "ringbuf.hpp"

using namespace std;

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
shared_ptr<T> Ringbuf<T>::push(T* t) {
    // Is there room?  (This tries to make room if not...)
    bool can = can_push();
    if (!can)
        return shared_ptr<T>();
    //cout << "Creating shared_ptr..." << endl;
    shared_ptr<T> p(t, _deleter);
    pthread_mutex_lock(&this->_live_lock);
    _live++;
    pthread_mutex_unlock(&this->_live_lock);
    pthread_mutex_lock(&this->_q_lock);
    _q.push_back(p);
    pthread_mutex_unlock(&this->_q_lock);
    //cout << "Now " << _live << " objects are live" << endl;
    return p;
}

template <class T>
shared_ptr<T> Ringbuf<T>::pop() {
    //cout << "Popping..." << endl;
    shared_ptr<T> p;
    pthread_mutex_lock(&this->_q_lock);
    if (_q.empty()) {
        //cout << "Pop: empty" << endl;
    } else {
        p = _q.pop_front();
        //cout << "Pop: returning " << *p << endl;
    }
    pthread_mutex_unlock(&this->_q_lock);
    return p;
}

template <class T>
vector<shared_ptr<T> > Ringbuf<T>::snapshot(bool (*testFunc)(const shared_ptr<T>)) {
    vector<shared_ptr<T> > vec;
    pthread_mutex_lock(&this->_q_lock);
    for (auto it = _q.begin(); it != _q.end(); it++) {
        if (!testFunc || testFunc(*it)) {
	    vec.push_back(*it);
        }
    }
    pthread_mutex_unlock(&this->_q_lock);
    return vec;
}

template <class T>
void Ringbuf<T>::dropping(shared_ptr<T> t) {}

// Called by the RingbufDeleter when a shared_ptr is deleted
template <class T>
void Ringbuf<T>::deleted(T* t) {
    //cout << "Deleting object: " << *t << endl;
    pthread_mutex_lock(&this->_live_lock);
    _live--;
    pthread_mutex_unlock(&this->_live_lock);
    //cout << "Now " << _live << " objects are live" << endl;
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
        //cout << "Push: _live >= _maxsize." << " (" << nlive << " >= " << _maxsize << ")" << endl;
        pthread_mutex_lock(&this->_q_lock);

        if (_q.empty()) {
            //cout << "Ring buffer empty but still too many live elements -- push fails" << endl;
            pthread_mutex_unlock(&this->_q_lock);
            return false;
        }
        //cout << "Dropping an element..." << endl;
        shared_ptr<T> p = _q.front();
        _q.pop_front();
        pthread_mutex_unlock(&this->_q_lock);
        dropping(p);
        p.reset();
        //cout << "Now " << _live << " live" << endl;
    }
    return true;
}

// Helper class that is the Deleter for our shared_ptr<frames>.  Calls
// Ringbuf.deleted() to track number of live frames.
template <class T>
RingbufDeleter<T>::RingbufDeleter(Ringbuf<T>* rb) : _ringbuf(rb) {}

template <class T>
void RingbufDeleter<T>::operator()(T* t) {
    //cout << "RingbufDelete::operator() called." << endl;
    _ringbuf->deleted(t);
}


