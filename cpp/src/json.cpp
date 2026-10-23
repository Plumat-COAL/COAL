#include "json.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace coal {

// ---------------------------------------------------------------- SHA-256 (FIPS 180-4)
namespace {
const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
inline uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
}  // namespace

Sha256::Sha256() {
    const uint32_t init[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                              0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::memcpy(h_, init, sizeof h_);
}

void Sha256::block(const unsigned char* p) {
    uint32_t w[64];
    for (int t = 0; t < 16; ++t)
        w[t] = (uint32_t)p[4 * t] << 24 | (uint32_t)p[4 * t + 1] << 16 | (uint32_t)p[4 * t + 2] << 8 | p[4 * t + 3];
    for (int t = 16; t < 64; ++t) {
        uint32_t s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >> 3);
        uint32_t s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >> 10);
        w[t] = w[t - 16] + s0 + w[t - 7] + s1;
    }
    uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3], e = h_[4], f = h_[5], g = h_[6], h = h_[7];
    for (int t = 0; t < 64; ++t) {
        uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = h + S1 + ch + K[t] + w[t];
        uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + mj;
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h_[0] += a; h_[1] += b; h_[2] += c; h_[3] += d; h_[4] += e; h_[5] += f; h_[6] += g; h_[7] += h;
}

void Sha256::update(const unsigned char* p, size_t n) {
    total_ += n;
    while (n) {
        size_t take = std::min(n, (size_t)64 - blen_);
        std::memcpy(buf_ + blen_, p, take);
        blen_ += take;
        p += take;
        n -= take;
        if (blen_ == 64) {
            block(buf_);
            blen_ = 0;
        }
    }
}

std::string Sha256::hex() {
    uint64_t bits = total_ * 8;
    unsigned char pad = 0x80;
    update(&pad, 1);
    unsigned char zero = 0;
    while (blen_ != 56) update(&zero, 1);
    unsigned char len[8];
    for (int k = 0; k < 8; ++k) len[k] = (unsigned char)(bits >> (56 - 8 * k));
    update(len, 8);
    char out[65];
    for (int k = 0; k < 8; ++k) std::snprintf(out + 8 * k, 9, "%08x", h_[k]);
    return std::string(out, 64);
}

std::string sha256_hex(const std::string& data) {
    Sha256 s;
    s.update((const unsigned char*)data.data(), data.size());
    return s.hex();
}

std::string sha256_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return "";
    Sha256 s;
    std::vector<char> buf(1 << 20);
    for (;;) {
        in.read(buf.data(), (std::streamsize)buf.size());
        std::streamsize n = in.gcount();
        if (n <= 0) break;
        s.update((const unsigned char*)buf.data(), (size_t)n);
    }
    return s.hex();
}

// ---------------------------------------------------------------- JSON
Json& Json::set(const std::string& k, Json v) {
    type = Object;
    for (auto& kv : o)
        if (kv.first == k) {
            kv.second = std::move(v);
            return kv.second;
        }
    o.push_back({k, std::move(v)});
    return o.back().second;
}

const Json* Json::get(const std::string& k) const {
    for (auto& kv : o)
        if (kv.first == k) return &kv.second;
    return nullptr;
}

static void dump_str(std::string& out, const std::string& s) {
    out += '"';
    for (unsigned char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            default:
                if (c < 0x20 || c >= 0x7f) {
                    char b[8];
                    std::snprintf(b, sizeof b, "\\u%04x", c);
                    out += b;
                } else {
                    out += char(c);
                }
        }
    }
    out += '"';
}

static void dump(std::string& out, const Json& j, int indent, int depth) {
    auto nl = [&](int d) {
        if (indent >= 0) {
            out += '\n';
            out.append((size_t)(indent * d), ' ');
        }
    };
    switch (j.type) {
        case Json::Null: out += "null"; break;
        case Json::Bool: out += j.b ? "true" : "false"; break;
        case Json::Int: out += std::to_string(j.i); break;
        case Json::Float: {
            char b[32];
            std::snprintf(b, sizeof b, "%.6g", j.f);
            out += b;
            if (!std::strpbrk(b, ".eEn")) out += ".0";
            break;
        }
        case Json::String: dump_str(out, j.s); break;
        case Json::Array:
            if (j.a.empty()) {
                out += "[]";
                break;
            }
            out += '[';
            for (size_t k = 0; k < j.a.size(); ++k) {
                if (k) out += indent >= 0 ? "," : ", ";
                nl(depth + 1);
                dump(out, j.a[k], indent, depth + 1);
            }
            nl(depth);
            out += ']';
            break;
        case Json::Object:
            if (j.o.empty()) {
                out += "{}";
                break;
            }
            out += '{';
            for (size_t k = 0; k < j.o.size(); ++k) {
                if (k) out += indent >= 0 ? "," : ", ";
                nl(depth + 1);
                dump_str(out, j.o[k].first);
                out += ": ";
                dump(out, j.o[k].second, indent, depth + 1);
            }
            nl(depth);
            out += '}';
            break;
    }
}

std::string json_dump(const Json& j, int indent) {
    std::string out;
    dump(out, j, indent, 0);
    return out;
}

void json_dump_to(std::string& out, const Json& j, int indent, int depth) { dump(out, j, indent, depth); }

namespace {
// Nesting deeper than this is malformed input (a checkpoint or report never nests more than a few levels):
// the parser recurses once per level, so without a limit a file of 100,000 '[' overflows the stack.
constexpr int MAX_DEPTH = 64;

struct Parser {
    const std::string& t;
    size_t p = 0;
    int depth = 0;
    [[noreturn]] void fail() { throw std::runtime_error("malformed JSON at offset " + std::to_string(p)); }
    void ws() {
        while (p < t.size() && (t[p] == ' ' || t[p] == '\n' || t[p] == '\r' || t[p] == '\t')) ++p;
    }
    bool lit(const char* w) {
        size_t n = std::strlen(w);
        if (t.compare(p, n, w) == 0) {
            p += n;
            return true;
        }
        return false;
    }
    std::string str() {
        if (t[p] != '"') fail();
        ++p;
        std::string s;
        while (p < t.size() && t[p] != '"') {
            char c = t[p++];
            if (c != '\\') {
                s += c;
                continue;
            }
            if (p >= t.size()) fail();
            char e = t[p++];
            switch (e) {
                case 'n': s += '\n'; break;
                case 'r': s += '\r'; break;
                case 't': s += '\t'; break;
                case 'b': s += '\b'; break;
                case 'f': s += '\f'; break;
                case 'u': {
                    if (p + 4 > t.size()) fail();
                    for (size_t k = p; k < p + 4; ++k)
                        if (!isxdigit((unsigned char)t[k])) fail();
                    unsigned v = (unsigned)std::stoul(t.substr(p, 4), nullptr, 16);
                    p += 4;
                    if (v < 0x80) s += char(v);
                    else if (v < 0x800) {
                        s += char(0xc0 | (v >> 6));
                        s += char(0x80 | (v & 63));
                    } else {
                        s += char(0xe0 | (v >> 12));
                        s += char(0x80 | ((v >> 6) & 63));
                        s += char(0x80 | (v & 63));
                    }
                    break;
                }
                default: s += e;
            }
        }
        if (p >= t.size()) fail();
        ++p;
        return s;
    }
    Json value() {
        if (++depth > MAX_DEPTH) fail();
        Json v = value_at_depth();
        --depth;
        return v;
    }
    Json value_at_depth() {
        ws();
        if (p >= t.size()) fail();
        char c = t[p];
        if (c == '{') {
            ++p;
            Json j = Json::object();
            ws();
            if (p < t.size() && t[p] == '}') {
                ++p;
                return j;
            }
            for (;;) {
                ws();
                std::string k = str();
                ws();
                if (p >= t.size() || t[p] != ':') fail();
                ++p;
                j.o.push_back({k, value()});
                ws();
                if (p < t.size() && t[p] == ',') {
                    ++p;
                    continue;
                }
                if (p < t.size() && t[p] == '}') {
                    ++p;
                    return j;
                }
                fail();
            }
        }
        if (c == '[') {
            ++p;
            Json j = Json::array();
            ws();
            if (p < t.size() && t[p] == ']') {
                ++p;
                return j;
            }
            for (;;) {
                j.a.push_back(value());
                ws();
                if (p < t.size() && t[p] == ',') {
                    ++p;
                    continue;
                }
                if (p < t.size() && t[p] == ']') {
                    ++p;
                    return j;
                }
                fail();
            }
        }
        if (c == '"') return Json(str());
        if (lit("true")) return Json(true);
        if (lit("false")) return Json(false);
        if (lit("null")) return Json();
        size_t q = p;
        bool flt = false;
        if (t[p] == '-') ++p;
        while (p < t.size() && (isdigit((unsigned char)t[p]) || t[p] == '.' || t[p] == 'e' || t[p] == 'E' ||
                                t[p] == '+' || t[p] == '-')) {
            if (!isdigit((unsigned char)t[p])) flt = true;
            ++p;
        }
        if (p == q) fail();
        std::string num = t.substr(q, p - q);
        try {
            if (flt) return Json(std::stod(num));
            return Json((long long)std::stoll(num));
        } catch (const std::exception&) {  // not a number, or out of range
            fail();
        }
    }
};
}  // namespace

Json json_parse(const std::string& text) {
    Parser ps{text};
    Json v = ps.value();
    ps.ws();
    if (ps.p != text.size()) ps.fail();
    return v;
}

}  // namespace coal
