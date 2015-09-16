/*
 * ICompressdFile.h
 *
 *  Created on: Apr 30, 2014
 *      Author: jrw32
 */

#ifndef UTILITY_ICOMPRESSEDFILE_H
#define UTILITY_ICOMPRESSEDFILE_H

#include <istream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/version.hpp>

// Workaround needed to access + change header_.state_ and header_.flags_
// NOTE: this is a REALLY bad idea!!
#if BOOST_VERSION < 104900

namespace private_access{
template <class Tag>
struct stowed
{
     static typename Tag::type value;
};
template <class Tag>
typename Tag::type stowed<Tag>::value;

// Generate a static data member whose constructor initializes
// stowed<Tag>::value.  This type will only be named in an explicit
// instantiation, where it is legal to pass the address of a private
// member.
template <class Tag, typename Tag::type x>
struct stow_private
{
     stow_private() { stowed<Tag>::value = x; }
     static stow_private instance;
};
template <class Tag, typename Tag::type x>
stow_private<Tag,x> stow_private<Tag,x>::instance;

}

struct hdr_flag {typedef int (boost::iostreams::detail::gzip_header::*type);};
struct hdr_state {typedef int (boost::iostreams::detail::gzip_header::*type);};

template class private_access::stow_private<hdr_flag, &boost::iostreams::detail::gzip_header::flags_>;
template class private_access::stow_private<hdr_state, &boost::iostreams::detail::gzip_header::state_>;
#endif


namespace boost{
namespace iostreams{



//------------------Definition of basic_bgzip_decompressor---------------------//

//
// Template name: basic_bgzip_decompressor
// Description: Model of InputFilter implementing compression in the
//      gzip format.
//
template<typename Alloc = std::allocator<char> >
class basic_bgzip_decompressor : basic_zlib_decompressor<Alloc> {
private:
    typedef basic_zlib_decompressor<Alloc>   base_type;
    typedef typename base_type::string_type  string_type;
public:
    typedef char char_type;
    struct category
        : dual_use,
          filter_tag,
          multichar_tag,
          closable_tag
        { };
    basic_bgzip_decompressor( int window_bits = gzip::default_window_bits,
                             int buffer_size = default_device_buffer_size );

    template<typename Sink>
    std::streamsize write(Sink& snk, const char_type* s, std::streamsize n)
    {
        std::streamsize result = 0;
        while(result < n) {
            if(state_ == s_start) {
                state_ = s_header;
                header_.reset();
                footer_.reset();
            }
            if (state_ == s_header) {
                int c = s[result++];
                header_.process(c);
                if (header_.done())
                    state_ = s_body;
            } else if (state_ == s_body) {
                try {
                    std::streamsize amt =
                        base_type::write(snk, s + result, n - result);
                    result += amt;
                    if (!this->eof()) {
                        break;
                    } else {
                        state_ = s_footer;
                    }
                } catch (const zlib_error& e) {
                    boost::throw_exception(gzip_error(e));
                }
            } else { // state_ == s_footer
                if (footer_.done()) {
                    if (footer_.crc() != this->crc())
                        boost::throw_exception(gzip_error(gzip::bad_crc));

                    base_type::close(snk, BOOST_IOS::out);
                    state_ = s_start;
                } else {
                    int c = s[result++];
                    footer_.process(c);
                }
            }
        }
        return result;
    }

    template<typename Source>
    std::streamsize read(Source& src, char_type* s, std::streamsize n)
    {
        typedef char_traits<char>  traits_type;
        std::streamsize            result = 0;
        peekable_source<Source>    peek(src, putback_);
        int hdr_bytes_read = 0;
        while (result < n && state_ != s_done) {
            if (state_ == s_start) {
                state_ = s_header;
                header_.reset();
                hdr_bytes_read = 0;
                footer_.reset();
            }
            if (state_ == s_header) {
                int c = boost::iostreams::get(peek);
                if (traits_type::is_eof(c)) {
                    boost::throw_exception(gzip_error(gzip::bad_header));
                } else if (traits_type::would_block(c)) {
                    break;
                }
                header_.process(c);
                ++hdr_bytes_read;
                // Workaround bug# 5908 for boost versions <= 1.48
#if BOOST_VERSION < 104900
				// after the 10th byte is read (the OS byte), we need to check
				// to see if the extra bit is set in header_.flags_
				// if it is, set the header_.state_ == s_xlen
				if(hdr_bytes_read == 10){

					if(header_.*private_access::stowed<hdr_flag>::value & gzip::flags::extra){
						--(header_.*private_access::stowed<hdr_state>::value);
					}
				}
#endif
                if (header_.done())
                    state_ = s_body;
            } else if (state_ == s_body) {
                try {
                    std::streamsize amt =
                        base_type::read(peek, s + result, n - result);
                    if (amt != -1) {
                        result += amt;
                        if (amt < n - result)
                            break;
                    } else {
                        peek.putback(this->unconsumed_input());
                        state_ = s_footer;
                    }
                } catch (const zlib_error& e) {
                    boost::throw_exception(gzip_error(e));
                }
            } else { // state_ == s_footer
                int c = boost::iostreams::get(peek);
                if (traits_type::is_eof(c)) {
                    boost::throw_exception(gzip_error(gzip::bad_footer));
                } else if (traits_type::would_block(c)) {
                    break;
                }
                footer_.process(c);
                if (footer_.done()) {
                    if (footer_.crc() != 0 && footer_.crc() != this->crc())
                        boost::throw_exception(gzip_error(gzip::bad_crc));
                    int c = boost::iostreams::get(peek);
                    if (traits_type::is_eof(c)) {
                        state_ = s_done;
                    } else {
                        peek.putback(c);
                        base_type::close(peek, BOOST_IOS::in);
                        state_ = s_start;
                        header_.reset();
                        footer_.reset();
                    }
                }
            }
        }
        if (peek.has_unconsumed_input()) {
            putback_ = peek.unconsumed_input();
        } else {
            putback_.clear();
        }
        return result != 0 || state_ != s_done ?
            result :
            -1;
    }

    template<typename Source>
    void close(Source& src, BOOST_IOS::openmode m)
    {
        try {
            base_type::close(src, m);
        } catch (const zlib_error& e) {
            state_ = s_start;
            boost::throw_exception(gzip_error(e));
        }
        if (m == BOOST_IOS::out) {
            if (state_ == s_start || state_ == s_header)
                boost::throw_exception(gzip_error(gzip::bad_header));
            else if (state_ == s_body)
                boost::throw_exception(gzip_error(gzip::bad_footer));
            else if (state_ == s_footer) {
                if (!footer_.done())
                    boost::throw_exception(gzip_error(gzip::bad_footer));
                else if(footer_.crc() != this->crc())
                    boost::throw_exception(gzip_error(gzip::bad_crc));
            } else {
                BOOST_ASSERT(!"Bad state");
            }
        }
        state_ = s_start;
    }

    std::string file_name() const { return header_.file_name(); }
    std::string comment() const { return header_.comment(); }
    bool text() const { return header_.text(); }
    int os() const { return header_.os(); }
    std::time_t mtime() const { return header_.mtime(); }
private:
    static gzip_params make_params(int window_bits);

    // Source adapter allowing an arbitrary character sequence to be put back.
    template<typename Source>
    struct peekable_source {
        typedef char char_type;
        struct category : source_tag, peekable_tag { };
        explicit peekable_source(Source& src, const string_type& putback = "")
            : src_(src), putback_(putback), offset_(0)
            { }
        std::streamsize read(char* s, std::streamsize n)
        {
            std::streamsize result = 0;

            // Copy characters from putback buffer
            std::streamsize pbsize =
                static_cast<std::streamsize>(putback_.size());
            if (offset_ < pbsize) {
                result = (std::min)(n, pbsize - offset_);
                BOOST_IOSTREAMS_CHAR_TRAITS(char)::copy(
                    s, putback_.data() + offset_, result);
                offset_ += result;
                if (result == n)
                    return result;
            }

            // Read characters from src_
            std::streamsize amt =
                boost::iostreams::read(src_, s + result, n - result);
            return amt != -1 ?
                result + amt :
                result ? result : -1;
        }
        bool putback(char c)
        {
            if (offset_) {
                putback_[--offset_] = c;
            } else {
                boost::throw_exception(
                    boost::iostreams::detail::bad_putback());
            }
            return true;
        }
        void putback(const string_type& s)
        {
            putback_.replace(0, offset_, s);
            offset_ = 0;
        }

        // Returns true if some characters have been putback but not re-read.
        bool has_unconsumed_input() const
        {
            return offset_ < static_cast<std::streamsize>(putback_.size());
        }

        // Returns the sequence of characters that have been put back but not re-read.
        string_type unconsumed_input() const
        {
            return string_type(putback_, offset_, putback_.size() - offset_);
        }
        Source&          src_;
        string_type      putback_;
        std::streamsize  offset_;
    };

    enum state_type {
        s_start   = 1,
        s_header  = s_start + 1,
        s_body    = s_header + 1,
        s_footer  = s_body + 1,
        s_done    = s_footer + 1
    };
    detail::gzip_header  header_;
    detail::gzip_footer  footer_;
    string_type          putback_;
    int                  state_;
};
BOOST_IOSTREAMS_PIPABLE(basic_bgzip_decompressor, 1)

typedef basic_bgzip_decompressor<> bgzip_decompressor;


//------------------Implementation of bgzip_decompressor-----------------------//

template<typename Alloc>
basic_bgzip_decompressor<Alloc>::basic_bgzip_decompressor
    (int window_bits, int buffer_size)
    : base_type(make_params(window_bits), buffer_size),
      state_(s_start)
    { }

template<typename Alloc>
gzip_params basic_bgzip_decompressor<Alloc>::make_params(int window_bits)
{
    gzip_params p;
    p.window_bits = window_bits;
    p.noheader = true;
    p.calculate_crc = true;
    return p;
}

//----------------------------------------------------------------------------//
}}

namespace BioBin {
namespace Utility {


/*!
 * \brief A class to read in compressed files automatically based on extension
 * This class, which has the same interface as an ifstream (so can be used
 * interchangeably) will automatically decompress
 */
class ICompressedFile : public std::istream {
public:
	ICompressedFile();
	explicit ICompressedFile(const char* fn, std::ios_base::openmode mode = std::ios_base::in);
	virtual ~ICompressedFile() {}

	void open(const char* fn, std::ios_base::openmode mode = std::ios_base::in);
	void close();

private:
	void openCompressed(const char* fn);

	boost::iostreams::filtering_istreambuf _infile;
	std::ifstream _base_f;


};

}

}

#endif /* ICOMPRESSDFILE_H_ */
