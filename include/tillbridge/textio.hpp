// TillBridge - reading whitespace- or comma-separated numeric text files.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace tb {

// Reads a text file line by line, skipping empty lines and lines starting with '#'.
class TextReader {
  public:
    explicit TextReader(const std::string& path) : path_(path), f_(std::fopen(path.c_str(), "rb")) {
        if (!f_) throw std::runtime_error("cannot open " + path);
    }
    ~TextReader() { if (f_) std::fclose(f_); }
    TextReader(const TextReader&) = delete;
    TextReader& operator=(const TextReader&) = delete;

    // Next non-empty, non-comment line; false at the end of the file.
    bool line(std::string& out) {
        out.clear();
        int c;
        for (;;) {
            out.clear();
            while ((c = std::fgetc(f_)) != EOF && c != '\n') out.push_back(char(c));
            if (c == EOF && out.empty()) return false;
            line_no_++;
            size_t a = out.find_first_not_of(" \t\r");
            if (a != std::string::npos && out[a] != '#') return true;
            if (c == EOF) return false;
        }
    }
    // Numbers of the next line (separators: blanks and commas). Fails at the end of the file.
    std::vector<double> numbers() {
        std::string s;
        if (!line(s)) fail("file ends early");
        return parse(s);
    }
    static std::vector<double> parse(const std::string& s) {
        std::vector<double> v;
        const char* p = s.c_str();
        for (;;) {
            while (*p == ' ' || *p == '\t' || *p == ',' || *p == '\r') p++;
            if (!*p) break;
            char* e;
            double d = std::strtod(p, &e);
            if (e == p) break;
            v.push_back(d);
            p = e;
        }
        return v;
    }
    [[noreturn]] void fail(const std::string& msg) const {
        throw std::runtime_error(msg + " in " + path_ + " (line " + std::to_string(line_no_) + ")");
    }

  private:
    std::string path_;
    std::FILE* f_;
    long line_no_ = 0;
};

}  // namespace tb
