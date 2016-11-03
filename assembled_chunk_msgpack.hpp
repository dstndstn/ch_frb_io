#ifndef _ASSEMBLED_CHUNK_MSGPACK_HPP
#define _ASSEMBLED_CHUNK_MSGPACK_HPP

#include <vector>
#include <iostream>

#include <msgpack.hpp>

#include <ch_frb_io.hpp>

/** Code for packing objects into msgpack mesages, and vice verse. **/

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
namespace adaptor {

// Place class template specialization here
template<>
struct convert<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    msgpack::object const& operator()(msgpack::object const& o,
                                      std::shared_ptr<ch_frb_io::assembled_chunk>& ch) const {
        if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
        //cout << "convert msgpack object to shared_ptr<assembled_chunk>..." << endl;
        if (o.via.array.size != 14) throw msgpack::type_error();
        msgpack::object* arr = o.via.array.ptr;

        // header string is [0]
        uint8_t version            = arr[1].as<uint8_t>();
        if (version != 1)
            throw std::runtime_error("ch_frb_io assembled_chunk msgpack version " + std::to_string(version) + ", expected 1");
        int beam_id                = arr[2].as<int>();
        int nupfreq                = arr[3].as<int>();
        int nt_per_packet          = arr[4].as<int>();
        int fpga_counts_per_sample = arr[5].as<int>();
        int nt_coarse              = arr[6].as<int>();
        int nscales                = arr[7].as<int>();
        int ndata                  = arr[8].as<int>();
        uint64_t ichunk            = arr[9].as<uint64_t>();
        uint64_t isample           = arr[10].as<uint64_t>();

        ch = ch_frb_io::assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
        ch->isample = isample;

        if (ch->nt_coarse != nt_coarse)
            throw std::runtime_error("ch_frb_l1 rpc: assembled_chunk nt_coarse mismatch");

        if (arr[11].type != msgpack::type::BIN) throw msgpack::type_error();
        if (arr[12].type != msgpack::type::BIN) throw msgpack::type_error();
        if (arr[13].type != msgpack::type::BIN) throw msgpack::type_error();

        int nsdata = nscales * sizeof(float);
        if (arr[11].via.bin.size != nsdata) throw msgpack::type_error();
        if (arr[12].via.bin.size != nsdata) throw msgpack::type_error();
        if (arr[13].via.bin.size != ndata ) throw msgpack::type_error();
        
        memcpy(ch->scales,  arr[11].via.bin.ptr, nsdata);
        memcpy(ch->offsets, arr[12].via.bin.ptr, nsdata);
        memcpy(ch->data,    arr[13].via.bin.ptr, ndata);

        return o;
    }
};

template<>
struct pack<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    template <typename Stream>
    packer<Stream>& operator()(msgpack::packer<Stream>& o, std::shared_ptr<ch_frb_io::assembled_chunk>  const& ch) const {
        // packing member variables as an array.
        //cout << "Pack shared_ptr<assembled-chunk> into msgpack object..." << endl;
        uint8_t version = 1;
        o.pack_array(14);
        o.pack("assembled_chunk in msgpack format");
        o.pack(version);
        o.pack(ch->beam_id);
        o.pack(ch->nupfreq);
        o.pack(ch->nt_per_packet);
        o.pack(ch->fpga_counts_per_sample);
        o.pack(ch->nt_coarse);
        o.pack(ch->nscales);
        o.pack(ch->ndata);
        o.pack(ch->ichunk);
        o.pack(ch->isample);
        // PACK FLOATS AS BINARY
        int nscalebytes = ch->nscales * sizeof(float);
        o.pack_bin(nscalebytes);
        o.pack_bin_body(reinterpret_cast<const char*>(ch->scales),
                        nscalebytes);
        o.pack_bin(nscalebytes);
        o.pack_bin_body(reinterpret_cast<const char*>(ch->offsets),
                        nscalebytes);

        o.pack_bin(ch->ndata);
        o.pack_bin_body(reinterpret_cast<const char*>(ch->data), ch->ndata);
        return o;
    }
};

    /* Apparently not needed yet?
template <>
struct object_with_zone<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    void operator()(msgpack::object::with_zone& o, std::shared_ptr<ch_frb_io::assembled_chunk>  const& v) const {
        o.type = type::ARRAY;
        cout << "Convert shared_ptr<assembled_chunk> into msgpack object_with_zone" << endl;
...
     */

} // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

#endif // _ASSEMBLED_CHUNK_MSGPACK_HPP
