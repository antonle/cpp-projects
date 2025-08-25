#include <iostream>
#include <cstring>

class String {
public:
    String(): str_(new char[1]), size_(0), capacity_(1) {}

    String(size_t counter, char symbol): str_(new char[counter]), size_(counter), capacity_(counter) {
        memset(str_, symbol, counter);
    }

    String(const char* other) {
        size_t size = strlen(other);
        str_ = new char[size];
        size_ = size;
        capacity_ = size;
        memcpy(str_, other, size_);
    }

    String(const String& other): String(other.size_, '\0') {
        memcpy(str_, other.str_, size_);
    }

    ~String() {
        clear();
    }

    size_t length() const {
        return size_;
    };

    void push_back(const char symbol) {
        if (str_ == nullptr) {
            remake();
        }

        if (size_ == capacity_) {
            increaseBuffer();
        }

        str_[size_++] = symbol;
    }

    void pop_back() {
        if (size_ * reduction_koef_ <= capacity_) {
            decreaseBuffer();
        }

        size_--;
    }

    char& front() {
        return str_[0];
    }

    const char& front() const {
        return str_[0];
    }

    char& back() {
        return str_[size_ - 1];
    }

    const char& back() const {
        return str_[size_ - 1];
    }

    size_t find(const String& other) const {
        int hasOverlap = 1;
        for (size_t i = 0; i < size_; ++i) {
            for (size_t j = 0; j < other.length(); ++j) {
                if (i + j >= size_ || str_[i + j] != other[j]) {
                    hasOverlap = 0;
                    break;
                }
            }
            if (hasOverlap == 1) {
                return i;
            }
            hasOverlap = 1;
        }
        return size_;
    }

    size_t rfind(const String& other) const {
        int hasOverlap = 1;
        for (size_t i = 0; i < size_; ++i) {
            for (size_t j = 0; j < other.length(); ++j) {
                if (size_ - 1 - i + j >= size_ || str_[size_ - 1 - i + j] != other[j]) {
                    hasOverlap = 0;
                    break;
                }
            }
            if (hasOverlap == 1) {
                return size_ - 1 - i;
            }
            hasOverlap = 1;
        }
        return size_;
    }

    String substr(size_t start, size_t count) const {
        String substring(count, '\0');
        for (size_t i = 0; i < count; ++i) {
            substring[i] = str_[start + i];
        }
        return substring;
    }

    bool empty() const {
        return !size_;
    }

    void clear() {
        delete[] str_;
        str_ = nullptr;
        size_ = 0;
        capacity_ = 1;
    }

    String& operator=(const String& other) {
        String copy = other;
        clear();

        size_ = copy.size_;
        capacity_ = copy.capacity_;
        str_ = new char[capacity_];
        memcpy(str_, copy.str_, size_);
        return *this;
    }

    char& operator[](size_t position) {
        return str_[position];
    }

    const char& operator[](size_t position) const {
        return str_[position];
    }

    String& operator+=(const String& rhs) {
        for (size_t i = 0; i < rhs.length(); ++i) {
            push_back(rhs[i]);
        }
        return *this;
    }

    String& operator+=(const char rhs) {
        push_back(rhs);
        return *this;
    }

private:
    void changeBuffer() {
        char* newstr = new char[capacity_];
        memcpy(newstr, str_, size_);
        delete[] str_;
        str_ = newstr;
    }

    void increaseBuffer() {
        capacity_ *= 2;
        changeBuffer();
    }

    void decreaseBuffer() {
        capacity_ /= 2;
        changeBuffer();
    }

    void remake() {
        capacity_ = 1;
        size_ = 0;
        str_ = new char[capacity_];
    }

    char* str_;
    size_t size_;
    size_t capacity_;
    static const int reduction_koef_ = 4;
};

bool operator==(const String& lhs, const String& rhs) {
    if (lhs.length() != rhs.length()) {
        return false;
    }
    for (size_t i = 0; i < lhs.length(); ++i) {
        if (lhs[i] != rhs[i]){
            return false;
        }
    }
    return true;
}

String operator+(char lhs, const String& rhs) {
    String copy(1, lhs);
    copy += rhs;
    return copy;
}

String operator+(const String& lhs, char rhs) {
    String copy = lhs;
    copy.push_back(rhs);
    return copy;
}

String operator+(const String& lhs, const String& rhs) {
    String copy = lhs;
    copy += rhs;
    return copy;
}

std::ostream& operator<<(std::ostream& out, const String& number) {
    for (size_t i = 0; i < number.length(); ++i){
        out << number[i];
    }
    return out;
}

std::istream& operator>>(std::istream& in, String& number) {
    number.clear();

    char symbol;
    while (in.get(symbol)) {
        if (symbol == '\n') {
            break;
        }
        number.push_back(symbol);
    }
    return in;
}

