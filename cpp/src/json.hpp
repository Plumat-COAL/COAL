// Minimal JSON (write + parse) and SHA-256, standard library only.
#pragma once
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace coal {

// ---------------------------------------------------------------- SHA-256
std::string sha256_hex(const std::string& data);
// Empty string if the file cannot be read.
std::string sha256_file(const std::string& path);

class Sha256 {
public:
    Sha256();
    void update(const unsigned char* p, size_t n);
    std::string hex();

private:
    void block(const unsigned char* p);
    uint32_t h_[8];
    unsigned char buf_[64];
    size_t blen_ = 0;
    uint64_t total_ = 0;
};

// ---------------------------------------------------------------- JSON values
struct Json {
    enum Type { Null, Bool, Int, Float, String, Array, Object } type = Null;
    bool b = false;
    long long i = 0;
    double f = 0;
    std::string s;
    std::vector<Json> a;
    std::vector<std::pair<std::string, Json>> o;  // insertion order kept

    Json() = default;
    static Json null() { return Json(); }
    Json(bool v) : type(Bool), b(v) {}
    Json(int v) : type(Int), i(v) {}
    Json(long long v) : type(Int), i(v) {}
    Json(size_t v) : type(Int), i((long long)v) {}
    Json(double v) : type(Float), f(v) {}
    Json(const char* v) : type(String), s(v) {}
    Json(const std::string& v) : type(String), s(v) {}
    static Json array() {
        Json j;
        j.type = Array;
        return j;
    }
    static Json object() {
        Json j;
        j.type = Object;
        return j;
    }
    Json& set(const std::string& k, Json v);  // object: add or replace
    Json& push(Json v) {
        a.push_back(std::move(v));
        return a.back();
    }
    const Json* get(const std::string& k) const;
    long long as_int(long long dflt = 0) const { return type == Int ? i : type == Float ? (long long)f : dflt; }
    std::string as_str() const { return type == String ? s : std::string(); }
};

// indent < 0: compact one line ("a": 1, separators like Python's json.dumps)
std::string json_dump(const Json& j, int indent = -1);
// Throws std::runtime_error on malformed input.
Json json_parse(const std::string& text);

}  // namespace coal
