#include <deque>

// forward decl
template<class T>
class RingbufDeleter;

template<class T>
class Ringbuf {
    friend class RingbufDeleter<T>;

public:
    /*
     Creates a new Ringbuf with the given maximum quota size.
     */
    Ringbuf(int maxsize);

    virtual ~Ringbuf() {}

    /*
     Try to push a new frame onto the ring buffer.  If the quota of
     frames has been reached, frames will be popped off the front
     until the buffer is empty.  If space is still not available, will
     return an empty shared_ptr.

     Returns the enqueued shared_ptr if successful.
     */
    shared_ptr<T> push(T* t);

    /*
     Check whether there is room to push a new frame onto the ring
     buffer.  This method will try to free up space by dropping frames
     off the front of the queue.  Returns false if no space can be found.
     */
    bool can_push();

    /*
     Pops the oldest frame from this buffer and returns it.  Returns
     an empty shared_ptr if the buffer is empty.
     */
    shared_ptr<T> pop();

    /*
     Retrieves frames from this ring buffer that satisfy the
     (optional) selection function.
     */
    vector<shared_ptr<T> > snapshot(bool (*testFunc)(const shared_ptr<T>));
    
protected:
    // Small helper class that calls deleted() when one of our frames
    // is deleted.
    RingbufDeleter<T> _deleter;

    // The actual queue of frames.
    deque<shared_ptr<T> > _q;

    // Number of live frames.
    size_t _live;
    // Maximum allowable number of live frames.
    size_t _maxsize;

    // We're about to drop this frame.
    virtual void dropping(shared_ptr<T> t) {}

    // Called by the RingbufDeleter when a shared_ptr is deleted
    void deleted(T* t);

};

// Helper class that is the Deleter for our shared_ptr<frames>.  Calls
// Ringbuf.deleted() to track number of live frames.
template<class T>
class RingbufDeleter {
public:
    RingbufDeleter(Ringbuf<T>* rb);
    void operator()(T* t);

protected:
    Ringbuf<T>* _ringbuf;
};

