
/**
 Sketch of the L1 ring buffer desired API -- what the RPC server would
 look like.

 ZMQ_DONTWAIT -- blocks on zmq_send when high-water-mark is reached (DEALER)
 */

//#include <zmq.hpp>


/*
 For write_chunks: we need to keep the assembled_chunks alive via a
 shared_ptr.  This can be achieved by the RPC server pulling requests
 off the client-facing socket in one thread, retrieving shared_ptrs to
 the assembled_chunks, and putting them in a (thread-safe) queue which
 will keep the shared_ptrs alive.  After enqueuing a chunk + args, it
 should send a request to a back-end worker to pull some work from the
 queue.  The worker can dequeue the assembled_chunk, msgpack it or
 otherwise write it to disk, and then drop the shared_ptr.
 */

/*
 For the telescoping time-sampling, we want to keep ~4 ring buffers.
 The property they have is that as the oldest data are dropped off the
 end of the ring buffer, we want to take two chunks and compile them
 into one chunk.  While this is happening, we still want to be able to
 access them. (ie, the RPC callbacks should still know about them).
 */

/*
 Each of these ~4 ring buffers should use a fixed amount of space; we
 may want to re-use the assembled_chunk memory.
 */

/*
 One idea: rather than explicitly "pinning" assembled_chunks, we could
 use shared_ptr.unique() -- if only the ring buffer has a copy of the
 shared_ptr, then that assembled_chunk can be dropped.
 */

/*
 Another idea is to use the assembled_chunk (or a wrapper of it)'s
 destructor to notify the ring buffer (or an allocator) that it is
 being deleted.
 */

/*
For the telescoping chain, we want the ring buffer to tell us when an
assembled_chunk is being dropped.
 */


/*
 assembled_chunk::make() is a factor for assembled_chunks.

 assembled_chunk_ringbuf::_make_assembled_chunk() wraps that or calls directly

 assembled_chunk_ringbuf holds active_chunk0 and active_chunk1
 *outside* the ringbuf.

 _put_assembled_chunk is where assembled_chunks get added to the ring buffer.

 */

/*
 assembled_chunk_ringbuf: needs to allocate a new chunk when it
 receives a packet with ichunk >= active_ichunk + 2.

 Could do:
 - try to allocate a new assembled_chunk.  If that fails:
 - if assembled_ringbuf_size == assembled_ringbuf_capacity,
   then the buffer is full of chunks waiting for the consumer thread.
   -> drop data?  This means we keep the oldest data; alternative
      would be to drop the oldest frame and keep the newest data
      instead.

 - iteratively: go to the next ringbuf position and drop that chunk
   off the back -- call callback, then reset the shared_ptr.  Then try
   to allocate a new assembled_chunk.  If resetting the shared_ptr
   drop the last reference to an assembled_chunk, then that chunk will
   become available.  If not (eg, if RPC thread has a copy of the
   shared_ptr), continue on to the next position.  This process will
   stop once we drop (assembled_ringbuf_capacity -
   assembled_ringbuf_size) chunks -- ie, all the chunks left in the
   buffer are the ones waiting to be consumed.

 - there will be a callback for when a chunk is being dropped.  This
   is where we handle the telescoping: we'll hold zero or one older
   chunks, and when a chunk is about to be dropped, if we had an older
   one, then we can downsample the two and then enqueue the resulting
   new chunk.  This enqueue call could (typically will, in
   steady-state) cause a cascade of dropping the oldest frame from the
   next ring buffer.


- note that the longer-duration ring buffers do not have a
  'downstream' consumer.

- Huh, what if the queue of assembled chunks waiting for the
  downstream consumer and the ring buffer per se are two separate data
  structures that each hold shared_ptrs?

 - ring buffer could allocate an assembled_chunk *subclass* that tells
   the ring buffer when it is deleted.

 */


/*
 Single ringbuffer public API.

 ringbuf<T>(int size)

 // called when an item is about to be dropped off the end.
 virtual void drop(T);

 // try to enqueue an element; returns false if the ring buffer is full.
 bool enqueue(const T&)
 bool enqueue(T*)?
 shared_ptr<T> enqueue(T*)?

 //
 vector<T> snapshot(bool test(T*)) ?




 where T : shared_ptr<assembled_chunk> and we will subclass it to
 override drop() to do the telescoping.
 

 */
#include <iostream>
#include <queue>

using namespace std;

// forward decl
template<class T>
class RingbufDeleter;

template<class T>
class Ringbuf {
    friend class RingbufDeleter<T>;

public:
    Ringbuf(int maxsize) :
        _deleter(this),
        _maxsize(maxsize)
    {}

    shared_ptr<T> push(T* t) {
        while (_live >= _maxsize) {
            cout << "Push: _live >= _maxsize." << endl;
            if (_q.empty()) {
                cout << "Ring buffer empty but still too many live elements -- push fails" << endl;
                //
                //return false;
                return shared_ptr<T>();
            }
            cout << "Dropping an element..." << endl;
            shared_ptr<T> p = _q.front();
            _q.pop();
            dropping(p);
        }
        cout << "Creating shared_ptr..." << endl;
        shared_ptr<T> p(t, _deleter);
        _live++;
        _q.push(p);
        //cout << "Pushing shared_ptr..." << endl;
        cout << "Now " << _live << " objects are live" << endl;
        return p;
    }

    shared_ptr<T> pop() {
        cout << "Popping..." << endl;
        shared_ptr<T> p = _q.front();
        _q.pop();
        cout << "Pop: returning " << *p << endl;
        return p;
    }
    
protected:
    RingbufDeleter<T> _deleter;

    // TEST
    std::queue<shared_ptr<T> > _q;

    size_t _live;
    size_t _maxsize;

    virtual void dropping(shared_ptr<T> t) {}

    // Called by the RingbufDeleter when a shared_ptr is deleted
    void deleted(T* t) {
        cout << "Deleting object: " << *t << endl;
        _live--;
        cout << "Now " << _live << " objects are live" << endl;
    }

};


template<class T>
class RingbufDeleter {
public:
    RingbufDeleter(Ringbuf<T>* rb) : _ringbuf(rb) {}

    void operator()(T* t) {
        cout << "RingbufDelete::operator() called." << endl;
        _ringbuf->deleted(t);
    }

protected:
    Ringbuf<T>* _ringbuf;
};

int main() {
    cout << "Creating ringbuf..." << endl;
    Ringbuf<int> rb(4);

    int a = 42;
    int b = 43;
    int c = 44;

    cout << "Pushing" << endl;
    rb.push(&a);
    cout << "Pushing" << endl;
    rb.push(&b);
    cout << "Pushing" << endl;
    rb.push(&c);

    cout << "Popping" << endl;
    shared_ptr<int> p1 = rb.pop();
    cout << "Popping" << endl;
    shared_ptr<int> p2 = rb.pop();
    cout << "Dropping" << endl;
    p1.reset();
    cout << endl;

    int d = 45;
    int e = 46;
    int f = 47;
    int g = 48;

    cout << "Pushing d..." << endl;
    shared_ptr<int> pd = rb.push(&d);

    cout << endl;
    cout << "Pushing e..." << endl;
    shared_ptr<int> pe = rb.push(&e);

    cout << endl;
    cout << "Pushing f..." << endl;
    shared_ptr<int> pf = rb.push(&f);

    cout << endl;
    cout << "Pushing g..." << endl;
    rb.push(&g);

    cout << "Done" << endl;

}


