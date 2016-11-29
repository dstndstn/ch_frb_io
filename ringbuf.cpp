
/**
 Sketch of the L1 ring buffer desired API -- what the RPC server would
 look like.

 */

/*
 For RPC write_chunks: we need to keep the assembled_chunks alive via
 a shared_ptr.  This can be achieved by the RPC server pulling
 requests off the client-facing socket in one thread, retrieving
 shared_ptrs to the assembled_chunks, and putting them in a
 (thread-safe) queue which will keep the shared_ptrs alive.  After
 enqueuing a chunk + args, it should send a request to a back-end
 worker to pull some work from the queue.  The worker can dequeue the
 assembled_chunk, msgpack it or otherwise write it to disk, and then
 drop the shared_ptr.

 ZMQ_DONTWAIT -- blocks on zmq_send when high-water-mark is reached (DEALER)
 //#include <zmq.hpp>
 */

/*
 For the telescoping time-sampling, we want to keep ~4 ring buffers.
 The property they have is that as the oldest data are dropped off the
 end of the ring buffer, we want to take two chunks and compile them
 into one chunk.  While this is happening, we still want to be able to
 access them. (ie, the RPC callbacks should still know about them).

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
 assembled_chunk::make() is a factory for assembled_chunks.

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


 - the longer-duration ring buffers do not have 'downstream'
   consumers.

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
#include <deque>

using namespace std;

#include "ringbuf.hpp"




template <class T>
Ringbuf<T>::Ringbuf(int maxsize) :
    _deleter(this),
    _q(),
    _live(0),
    _maxsize(maxsize)
{}

template <class T>
Ringbuf<T>::~Ringbuf() {}

template <class T>
shared_ptr<T> Ringbuf<T>::push(T* t) {
    bool can = _can_push();
    if (!can)
        return shared_ptr<T>();


    cout << "Creating shared_ptr..." << endl;
    shared_ptr<T> p(t, _deleter);
    _live++;
    _q.push_back(p);
    cout << "Now " << _live << " objects are live" << endl;
    return p;
}

template <class T>
bool Ringbuf<T>::can_push() {
    return _can_push();
}

template <class T>
shared_ptr<T> Ringbuf<T>::pop() {
    cout << "Popping..." << endl;
    if (_q.empty()) {
        cout << "Pop: empty" << endl;
        return shared_ptr<T>();
    }
    shared_ptr<T> p = _q.pop_front();
    cout << "Pop: returning " << *p << endl;
    return p;
}

template <class T>
vector<shared_ptr<T> > Ringbuf<T>::snapshot(bool (*testFunc)(const shared_ptr<T>)) {
    vector<shared_ptr<T> > vec;
    for (auto it = _q.begin(); it != _q.end(); it++) {
        if (!testFunc || testFunc(*it)) {
	    vec.push_back(*it);
        }
    }
    return vec;
}

template <class T>
void Ringbuf<T>::dropping(shared_ptr<T> t) {}

// Called by the RingbufDeleter when a shared_ptr is deleted
template <class T>
void Ringbuf<T>::deleted(T* t) {
    //cout << "Deleting object: " << *t << endl;
    _live--;
    cout << "Now " << _live << " objects are live" << endl;
    // FIXME --?
    delete t;
}

template <class T>
bool Ringbuf<T>::_can_push() {
    while (_live >= _maxsize) {
        cout << "Push: _live >= _maxsize." << " (" << _live << " >= " << _maxsize << ")" << endl;
        if (_q.empty()) {
            cout << "Ring buffer empty but still too many live elements -- push fails" << endl;
            return false;
        }
        cout << "Dropping an element..." << endl;
        shared_ptr<T> p = _q.front();
        _q.pop_front();
        dropping(p);
        p.reset();
        cout << "Now " << _live << " live" << endl;
    }
    return true;
}

// Helper class that is the Deleter for our shared_ptr<frames>.  Calls
// Ringbuf.deleted() to track number of live frames.
template <class T>
RingbufDeleter<T>::RingbufDeleter(Ringbuf<T>* rb) : _ringbuf(rb) {}

template <class T>
void RingbufDeleter<T>::operator()(T* t) {
    cout << "RingbufDelete::operator() called." << endl;
    _ringbuf->deleted(t);
}



////////////////////////// How we'd use this ring buffer in L1 //////////////////////////


#include "ch_frb_io.hpp"
using namespace ch_frb_io;

std::ostream& operator<<(std::ostream& s, const assembled_chunk& ch) {
    s << "assembled_chunk(beam " << ch.beam_id << ", ichunk " << ch.ichunk << ")";
    return s;
}

class L1Ringbuf;

class AssembledChunkRingbuf : public Ringbuf<assembled_chunk> {

public:
    AssembledChunkRingbuf(int binlevel, L1Ringbuf* parent, int maxsize) :
        Ringbuf<assembled_chunk>(maxsize),
        _binlevel(binlevel),
        _parent(parent)
    {}

    virtual ~AssembledChunkRingbuf() {}

protected:
    // my time-binning. level: 0 = original intensity stream; 1 =
    // binned x 2, 2 = binned x 4.
    int _binlevel;
    L1Ringbuf* _parent;

    // Called when the given frame *t* is being dropped off the buffer
    // to free up some space for a new frame.
    virtual void dropping(shared_ptr<assembled_chunk> t);

};

class L1Ringbuf {
    friend class AssembledChunkRingbuf;

    static const size_t Nbins = 4;

public:
    L1Ringbuf() :
        _q(),
        _rb(),
        _dropped()
    {
        // Create the ring buffer objects for each time binning
        // (0 = native rate, 1 = binned by 2, ...)
        for (size_t i=0; i<Nbins; i++)
            _rb.push_back(shared_ptr<AssembledChunkRingbuf>
                          (new AssembledChunkRingbuf(i, this, 4)));
        // Fill the "_dropped" array with empty shared_ptrs.
        for (size_t i=0; i<Nbins-1; i++)
            _dropped.push_back(shared_ptr<assembled_chunk>());
    }

    /*
     Tries to enqueue an assembled_chunk.  If no space can be
     allocated, returns false.  The ring buffer now assumes ownership
     of the assembled_chunk.
     */
    bool push(assembled_chunk* ch) {
        shared_ptr<assembled_chunk> p = _rb[0]->push(ch);
        if (!p)
            return false;
        _q.push_back(p);
        return true;
    }

    /*
     Returns the next assembled_chunk for downstream processing.
     */
    shared_ptr<assembled_chunk> pop() {
        if (_q.empty())
            return shared_ptr<assembled_chunk>();
        shared_ptr<assembled_chunk> p = _q.front();
        _q.pop_front();
        return p;
    }

    /*
     Prints a report of the assembled_chunks currently queued.
     */
    void print() {
        cout << "L1 ringbuf:" << endl;
        cout << "  downstream: [ ";
        for (auto it = _q.begin(); it != _q.end(); it++) {
            cout << (*it)->ichunk << " ";
        }
        cout << "];" << endl;
        for (size_t i=0; i<Nbins; i++) {
            vector<shared_ptr<assembled_chunk> > v = _rb[i]->snapshot(NULL);
            cout << "  binning " << i << ": [ ";
            for (auto it = v.begin(); it != v.end(); it++) {
                cout << (*it)->ichunk << " ";
            }
            cout << "]" << endl;
            if (i < Nbins-1) {
                cout << "  dropped " << i << ": ";
                if (_dropped[i])
                    cout << _dropped[i]->ichunk << endl;
                else
                    cout << "none" << endl;
            }
        }
    }
    
protected:
    // The queue for downstream
    deque<shared_ptr<assembled_chunk> > _q;

    // The ring buffers for each time-binning.  Length fixed at Nbins.
    vector<shared_ptr<AssembledChunkRingbuf> > _rb;

    // The assembled_chunks that have been dropped from the ring
    // buffers and are waiting for a pair to be time-downsampled.
    // Length fixed at Nbins-1.
    vector<shared_ptr<assembled_chunk> > _dropped;

    // Called from the AssembledChunkRingbuf objects when a chunk is
    // about to be dropped from one binning level of the ringbuf.  If
    // the chunk does not have a partner waiting (in _dropped), then
    // it is saved in _dropped.  Otherwise, the two chunks are merged
    // into one new chunk and added to the next binning level's
    // ringbuf.
    void dropping(int binlevel, shared_ptr<assembled_chunk> ch) {
        cout << "Bin level " << binlevel << " dropping a chunk" << endl;
        if (binlevel >= (int)(Nbins-1))
            return;

        if (_dropped[binlevel]) {
            cout << "Now have 2 dropped chunks from bin level " << binlevel << endl;
            // FIXME -- bin down
            assembled_chunk* binned = new assembled_chunk(ch->beam_id, ch->nupfreq, ch->nt_per_packet, ch->fpga_counts_per_sample, _dropped[binlevel]->ichunk);
            // push onto _rb[level+1]
            _rb[binlevel+1]->push(binned);
            _dropped[binlevel].reset();
        } else {
            // Keep this one until its partner arrives!
            cout << "Saving as _dropped" << binlevel << endl;
            _dropped[binlevel] = ch;
        }
    }

};

// after L1Ringbuf has been declared...
void AssembledChunkRingbuf::dropping(shared_ptr<assembled_chunk> t) {
    _parent->dropping(_binlevel, t);
}

#include <new>
#include <memory>

class MyClass {
public:
    MyClass(int x) : _x(x) {
        cout << "MyClass(" << x << ")" << endl;
    }
    ~MyClass() {
        cout << "~MyClass(" << _x << ")" << endl;
    }

    static int _nlive;

    static void* operator new(size_t count) {
        _nlive++;
        cout << "MyClass.new(" << count << "); nlive=" << _nlive << endl;
        return ::operator new(count);
    }

    static void operator delete(void* ptr) {
        _nlive--;
        cout << "MyClass.delete(); nlive=" << _nlive << endl;
        ::operator delete(ptr);
    }

    int _x;

};

int MyClass::_nlive = 0;

int main() {

    cout << "Create" << endl;

    MyClass a(42);

    MyClass* b = new MyClass(43);

    // Interestingly, this does NOT call the MyClass::operator new
    // (because make_shared merges the shared_ptr control block and
    // MyClass memory)
    shared_ptr<MyClass> c = make_shared<MyClass>(44);

    // make_unique is new in C++14
    //unique_ptr<MyClass> d = std::make_unique<MyClass>(45);
    unique_ptr<MyClass> d(new MyClass(45));

    cout << "Delete" << endl;

    d.reset();
    c.reset();
    delete b;

    exit(0);


    L1Ringbuf rb;

    int beam = 77;
    int nupfreq = 4;
    int nt_per = 16;
    int fpga_per = 400;

    assembled_chunk* ch;
    //ch = assembled_chunk::make(4, nupfreq, nt_per, fpga_per, 42);

    std::random_device rd;
    std::mt19937 rng(rd());
    rng.seed(42);
    std::uniform_int_distribution<> rando(0,1);

    for (int i=0; i<100; i++) {
        ch = new assembled_chunk(beam, nupfreq, nt_per, fpga_per, i);
        rb.push(ch);

        cout << "Pushed " << i << endl;
        rb.print();
        cout << endl;

        // downstream thread consumes with a lag of 2...
        if (i >= 2) {
            // Randomly consume 0 to 2 chunks
            if (rando(rng)) {
                cout << "Downstream consumes a chunk" << endl;
                rb.pop();
            }
            if (rando(rng)) {
                cout << "Downstream consumes a chunk" << endl;
                rb.pop();
            }
        }
    }

}


/*
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
 */
