//Created by Yury Vazhenin
#include <iostream>
#include <compare>
#include <vector>

template <size_t N, size_t D>
struct is_prime_rec {
    static constexpr bool value = (N % D == 0) ? false : (is_prime_rec<N, D - 1>::value);
};

template<size_t N>
struct is_prime_rec<N, 1> {
    static constexpr bool value = true;
};

template<size_t N>
struct is_prime_rec<N, 0> {
    static constexpr bool value = false;
};

constexpr int get_sqrt_rec(const size_t N, const size_t D) {
    if (D * D < N) {
        return get_sqrt_rec(N, D + 1);
    }
    return D;
}

template<size_t N>
struct is_prime {
    static constexpr bool value = is_prime_rec<N, get_sqrt_rec(N, 0)>::value;
};

template<>
struct is_prime<1> {
    static constexpr bool value = false;
};


long long bin_pow (long long a, int n, long long mod) {
    if (n == 0)
        return 1;
    if (n % 2) {
        return (a * bin_pow(a, n - 1, mod)) % mod;
    } else {
        long long ans = bin_pow(a, n / 2, mod);
        return (ans * ans) % mod;
    }
}

template<int N>
class Residue {
    int n;
public:
    Residue() : n(0) {}

    Residue(long long n) : n(((n % N) + N) % N) {}

    Residue& operator*=(const Residue& r) {
        n = (n * r.n) % N;
        return *this;
    }

    Residue operator*(const Residue& r) const {
        return Residue(n * r.n);
    }

    Residue& operator+=(const Residue& r) {
        n = (n + r.n) % N;
        return *this;
    }

    Residue operator+(const Residue& r) const {
        return Residue(n + r.n);
    }

    Residue& operator-=(const Residue& r) {
        n = (N + n - r.n) % N;
        return *this;
    }

    Residue operator-(const Residue& r) const {
        return Residue(n - r.n + N);
    }

    bool operator==(const Residue& r) const = default;

    template<int M = N>
    std::enable_if_t<is_prime<M>::value, Residue&> operator/=(const Residue& r) {
        n = (((long long) n * bin_pow(r.n, N - 2, N)) % (long long) N);
        return *this;
    }

    explicit operator int() {
        return n;
    }

    ~Residue() = default;
};



template<int N>
std::enable_if_t<is_prime<N>::value, Residue<N>> operator/ (const Residue<N>& a, const Residue<N>& b) {
    Residue<N> answer = a;
    answer /= b;
    return answer;
}


//BigInteger start

enum Sign {
    plus,
    minus
};

class BigInteger {
private:
    int size, capacity;
    Sign sign = plus;
    char* arr;
public:
    BigInteger();

    BigInteger(long long x);

    BigInteger(const BigInteger& number);

    BigInteger(int size, int capacity, Sign sign, char* arr);

    ~BigInteger();

    BigInteger& operator=(const BigInteger& other);

    void multiply(int deg);

    Sign get_sign() const;

    void change_sign();

    BigInteger abs() const;

    int get_size() const;

    bool is_zero() const;

    void push_back(char next_symbol);

    void reverse();

    void clear();

    BigInteger& operator*=(const BigInteger& other);

    BigInteger& operator/=(const BigInteger& other);

    BigInteger& operator+=(const BigInteger& other);

    BigInteger& operator-=(const BigInteger& other);

    BigInteger& operator%=(const BigInteger& n);

    BigInteger operator-() const;

    std::string toString() const;

    explicit operator bool() const;

    BigInteger& operator++();

    BigInteger operator++(int);

    BigInteger& operator--();

    BigInteger operator--(int);

    explicit operator double() const;

    friend std::strong_ordering operator<=>(const BigInteger&, const BigInteger&);
private:
    BigInteger(Sign _sign, int _capacity); //не могу удалить этот конструктор, все же он довольно часто полезен
    // , чтобы выделить сразу необходимое количество памяти. тем не менее, сделал его приватным, чтобы у юзера ене было к нему доступа

    BigInteger& minus_one();

    BigInteger& plus_one();

    void increase_capacity(int new_size);

    void cut();

    const char* data() const;

    char* data();

    static void karatsuba(int*, int*, int*, int);

    static void initialize_arrays(BigInteger&, const BigInteger&, const char*&, const char*&);

    static std::pair<int*, int*> get_two_arrays(const BigInteger&, const BigInteger&);

    int quotient(BigInteger*, BigInteger&);

    char* finalize(int*, int);

    void addition(const BigInteger& other);

    void subtraction(const BigInteger& other);
};


BigInteger::BigInteger() : size(1), capacity(1), arr(new char[1] {'\0'}) {}

BigInteger::BigInteger(long long number) {
    if (number < 0ll) {
        sign = minus;
        number = -number;
    }
    if (number == 0ll) {
        size = 2; capacity = 2;
        arr = new char[capacity];
        arr[0] = '0';
        arr[1] = '\0';
        return;
    }
    int current = 0;
    capacity = 20;
    arr = new char[capacity];
    while (number != 0) {
        arr[current++] = char('0' + number % 10);
        number /= 10;
    }
    arr[current++] = '\0';
    size = current;
}

BigInteger::BigInteger(Sign _sign, int _capacity) : size(1), capacity(_capacity), sign(_sign), arr(new char[capacity]){}

BigInteger::BigInteger(int size, int capacity, Sign sign, char* arr) : size(size), capacity(capacity), sign(sign), arr(arr) {}

BigInteger::BigInteger(const BigInteger& other) : size(other.size), capacity(other.capacity)
        , sign(other.sign), arr(new char[capacity]) {
    memcpy(arr, other.arr, size * sizeof(char));
}

BigInteger::~BigInteger() {
    delete[] arr;
}

BigInteger& BigInteger::operator=(const BigInteger& other) {
    if (this == &other)
        return *this;
    sign = other.sign;
    if (capacity > other.size) {
        size = other.size;
        memcpy(arr, other.arr, size * sizeof(char));
        return *this;
    }
    size = other.size;
    capacity = std::max(other.size, capacity * 2);
    delete[] arr;
    arr = new char[capacity];
    memcpy(arr, other.arr, size * sizeof(char));
    return *this;
}

Sign BigInteger::get_sign() const {
    return sign;
}

const char* BigInteger::data() const {
    return arr;
}

char* BigInteger::data() {
    return arr;
}

void BigInteger::change_sign() {
    sign = Sign(1 - sign);
}

int BigInteger::get_size() const {
    return size - 1;
}

bool BigInteger::is_zero() const {
    return arr[get_size() - 1] == '0';
}

std::strong_ordering operator<=>(const BigInteger& a, const BigInteger& b) {
    if (b.is_zero() && a.is_zero())
        return std::strong_ordering::equal;
    if (a.is_zero())
        return (b.get_sign() == plus ? std::strong_ordering::less : std::strong_ordering::greater);
    if (b.is_zero())
        return (a.get_sign() == plus ? std::strong_ordering::greater : std::strong_ordering::less);
    if (a.get_sign() != b.get_sign())
        return (a.get_sign() == plus ? std::strong_ordering::greater : std::strong_ordering::less);
    if (a.get_size() != b.get_size())
        return ((a.get_size() < b.get_size()) ^ (a.get_sign() == plus) ? std::strong_ordering::greater : std::strong_ordering::less);
    for (int i = a.get_size() - 1; i >= 0; i--) {
        if (a.data()[i] != b.data()[i])
            return ((a.get_sign() == plus) ^ (a.data()[i] < b.data()[i]) ? std::strong_ordering::greater : std::strong_ordering::less);
    }
    return std::strong_ordering::equal;
}

bool operator==(const BigInteger& a, const BigInteger& b) {
    return (a <=> b) == std::strong_ordering::equal;
}
void BigInteger::increase_capacity(int new_size) {
    capacity = std::max(capacity * 2, new_size);
    char* new_arr = new char[capacity];
    memcpy(new_arr, arr, size * sizeof(char));
    delete[] arr;
    arr = new_arr;
}

void BigInteger::addition(const BigInteger &other) {
    int first_size = get_size(), second_size = other.get_size();
    if (std::max(first_size, second_size) + 2 > capacity)
        increase_capacity(std::max(first_size, second_size) + 2);

    int temporary = 0;
    size = std::max(first_size, second_size);
    for (int i = 0; i < size; i++) {
        if (i < first_size)
            temporary += (arr[i] - '0');
        if (i < second_size)
            temporary += (other.arr[i] - '0');
        arr[i] = char('0' + temporary % 10);
        temporary /= 10;
    }
    if (temporary)
        arr[size++] = char('0' + temporary);
    arr[size++] = '\0';
}

BigInteger& BigInteger::operator+=(const BigInteger& other) {
    if (other == 0)
        return *this;
    if (other == 1)
        return ++*this;
    if (sign != other.sign) {
        change_sign();
        subtraction(other);
        change_sign();
    } else {
        addition(other);
    }
    return *this;
}

void BigInteger::initialize_arrays(BigInteger& a, const BigInteger& b, const char*& arr1, const char*& arr2) {
    if (a.get_sign() == minus) {
        if (a < b) {
            arr1 = a.data();
            arr2 = b.data();
        } else {
            a.change_sign();
            arr1 = b.data();
            arr2 = a.data();
        }
    } else {
        if (a < b) {
            a.change_sign();
            arr1 = b.data();
            arr2 = a.data();
        } else {
            arr1 = a.data();
            arr2 = b.data();
        }
    }
}

void BigInteger::subtraction(const BigInteger &other) {
    const char* arr1;
    const char* arr2;
    BigInteger::initialize_arrays(*this, other, arr1, arr2);

    int first_size = std::max(get_size(), other.get_size());
    int second_size = std::min(get_size(), other.get_size());
    bool mark[first_size];
    std::fill(mark, mark + first_size, false);
    size = first_size;
    for (int i = 0; i < first_size; i++) {
        int temporary = 0;
        temporary += (10 + (arr1[i] - '0') - int(mark[i])) % 10;
        if (i >= second_size) {
            arr[i] = char('0' + temporary);
            continue;
        } else if (temporary >= (arr2[i] - '0')) {
            arr[i] = char('0' + temporary - (arr2[i] - '0'));
            continue;
        }
        int in = i + 1;
        while (arr1[in] == '0') {
            mark[in] = true;
            in++;
        }
        mark[in] = true;
        temporary += 10;
        arr[i] = char('0' + temporary - (arr2[i] - '0'));
    }
    arr[size++] = '\0';
    cut();
}

BigInteger& BigInteger::operator-=(const BigInteger& other) {
    if (other == 0)
        return *this;
    if (other == 1)
        return --*this;
    if (sign != other.sign) {
        change_sign();
        addition(other);
        change_sign();
    } else {
        if (capacity < std::max(size, other.size))
            increase_capacity(std::max(size, other.size));
        subtraction(other);
    }
    return *this;
}

BigInteger operator+(const BigInteger& n, const BigInteger& other) {
    BigInteger answer = n;
    answer += other;
    return answer;
}

BigInteger operator-(const BigInteger& n, const BigInteger& other) {
    BigInteger answer = n;
    answer -= other;
    return answer;
}

void BigInteger::cut() {
    int in = size - 2;
    while (in > 0 && arr[in] == '0') {
        arr[in] = '\0';
        in--;
        size--;
    }
}

void BigInteger::multiply(int deg) {
    if (!*this)
        return;
    if (size + deg > capacity) {
        capacity = size + deg + 2;
        char *new_arr = new char[capacity];
        memcpy(new_arr + deg, arr, size * sizeof(char));
        std::fill(new_arr, new_arr + deg, '0');
        delete[] arr;
        arr = new_arr;
        size = size + deg;
    } else {
        for (int i = size - 1; i >= 0; i--) {
            arr[i + deg] = arr[i];
        }
        std::fill(arr, arr + deg, '0');
        size += deg;
    }
}

void BigInteger::karatsuba(int* first_number, int* second_number, int* result, int number_size) {
    if (number_size <= 128) {
        for (int i = 0; i < number_size; i++) {
            for (int j = 0; j < number_size; j++) {
                result[i + j] += (first_number[i] * second_number[j]);
            }
        }
        return;
    }
    int mid = number_size / 2;
    int sum1[mid], sum2[mid], product_of_two_sums[number_size];
    std::fill(product_of_two_sums, product_of_two_sums + number_size, 0);
    for (int i = 0; i < mid; i++) {
        sum1[i] = first_number[i] + first_number[i + mid];
        sum2[i] = second_number[i] + second_number[i + mid];
    }
    BigInteger::karatsuba(sum1, sum2, product_of_two_sums, mid);
    BigInteger::karatsuba(first_number, second_number, result, mid);
    BigInteger::karatsuba(first_number + mid, second_number + mid, result + number_size, mid);
    for (int i = 0; i < mid; i++) {
        int c1 = result[i + mid] + product_of_two_sums[i] - result[i] - result[i + number_size];
        int c2 = result[i + number_size] + product_of_two_sums[i + mid] - result[i + mid] - result[i + number_size + mid];
        result[i + mid] = c1, result[i + number_size] = c2;
    }
}

std::pair<int*, int*> BigInteger::get_two_arrays (const BigInteger& a, const BigInteger& b) {
    int size1 = 1 << int(ceil(log2(a.get_size())));
    int size2 = 1 << int(ceil(log2(b.get_size())));
    int* first_number = new int[std::max(size1, size2)] {};
    int* second_number = new int[std::max(size1, size2)] {};
    for (int i = 0; i < a.get_size(); i++)
        first_number[i] = a.data()[i] - '0';
    for (int i = 0; i < b.get_size(); i++)
        second_number[i] = b.data()[i] - '0';
    return std::make_pair(first_number, second_number);
}

char* BigInteger::finalize(int* result, int sz) {
    char* final_result = new char[2 * sz + 1];
    std::fill(final_result, final_result + 2 * sz, '\0');
    for (int i = 0; i < 2 * sz + 1; i++) {
        if (i != 2 * sz) {
            result[i + 1] += result[i] / 10;
        }
        final_result[i] = char('0'+ (result[i] % 10));
    }
    return final_result;
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
    int sz = 1 << int(ceil(log2(std::max(get_size(), other.get_size()))));
    std::pair<int*, int*> pair_of_arrays = BigInteger::get_two_arrays(*this, other);
    int* result = new int[2 * sz + 1];
    std::fill(result, result + 2 * sz, 0);
    BigInteger::karatsuba(pair_of_arrays.first, pair_of_arrays.second, result, sz);
    char* final_result = finalize(result, sz);
    delete[] arr;
    arr = final_result;
    size = capacity = 2 * sz + 1;
    if (other.sign != sign) {
        sign = minus;
    } else {
        sign = plus;
    }
    cut();
    delete[] pair_of_arrays.first;
    delete[] pair_of_arrays.second;
    delete[] result;
    return *this;
}

int BigInteger::quotient (BigInteger* help, BigInteger& now) {
    int l = 0, r = 10;
    while (l + 1 < r) {
        int cur = (l + r) / 2;
        if (help[cur] <= now) {
            l = cur;
        } else {
            r = cur;
        }
    }
    return l;
}

BigInteger& BigInteger::operator/=(const BigInteger& other) {
    BigInteger other_abs = other.abs();
    Sign final_sign = (sign == other.sign ? plus : minus);
    sign = plus;
    if (*this < other_abs) {
        *this = 0;
        return *this;
    }
    BigInteger answer(final_sign, 2);
    BigInteger now (plus, other.size +  1);
    memcpy(now.arr, arr + size - other.size, other.size * sizeof(char));
    int left = size - other.size;
    now.size = other.size;
    if (now < other_abs) {
        now.multiply(1);
        now.arr[0] = arr[--left];
    }
    BigInteger help[10];
    help[0] = 0;
    for (int i = 1; i < 10; i++) {
        help[i] = help[i - 1] + other_abs;
    }
    while (true) {
        int quot = quotient(help, now);
        answer.push_back(char('0' + quot));
        now -= (help[quot]);
        if (left == 0) {
            break;
        }
        now.multiply(1);
        now.arr[0] = arr[--left];
        now.cut();
    }
    answer.reverse();
    *this = answer;
    return *this;
}

BigInteger BigInteger::abs() const{
    BigInteger answer = *this;
    answer.sign = plus;
    return answer;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
    BigInteger answer = a;
    answer *= b;
    return answer;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b) {
    BigInteger answer = a;
    answer /= b;
    return answer;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
    *this -= (*this / other) * other;
    return *this;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
    BigInteger answer = a - (a / b) * b;
    return answer;
}

BigInteger BigInteger::operator-() const{
    BigInteger result = *this;
    result.change_sign();
    return result;
}

BigInteger::operator bool() const {
    return arr[get_size() - 1] != '0';
}

std::string BigInteger::toString() const {
    if (arr[get_size() - 1] == '0') {
        char answer[2];
        answer[0] = '0';
        answer[1] = '\0';
        return answer;
    }
    if (sign == plus) {
        char answer[size];
        for (int i = 0; i < size - 1; i++) {
            answer[i] = arr[get_size() - 1 - i];
        }
        answer[size - 1] = '\0';
        return answer;
    }
    char answer[size + 1];
    answer[0] = '-';
    for (int i = 0; i < size - 1; i++) {
        answer[1 + i] = arr[get_size() - 1 - i];
    }
    answer[size] = '\0';
    return answer;
}

void BigInteger::clear() {
    sign = plus;
    arr[0] = '\0';
    size = 1;
}

void BigInteger::push_back(char next_symbol) {
    if (size == capacity) {
        capacity *= 2;
        char* new_arr = new char[capacity];
        memcpy(new_arr, arr, size * sizeof(char));
        delete[] arr;
        arr = new_arr;
    }
    arr[size - 1] = next_symbol;
    arr[size++] = '\0';
}

void BigInteger::reverse() {
    for (int i = 0; i < (get_size() / 2); i++) {
        std::swap(arr[i], arr[size - 2 - i]);
    }
}

BigInteger& BigInteger::minus_one() {
    int now = 0;
    while ((now < size) && (arr[now] == '0')) {
        arr[now++] = '9';
    }
    arr[now]--;
    cut();
    return *this;
}

BigInteger& BigInteger::plus_one() {
    int now = 0;
    while ((now < size) && (arr[now] == '9')) {
        arr[now] = '0';
        now++;
    }
    if (now == size - 1) {
        push_back('1');
    } else {
        arr[now]++;
    }
    return *this;
}

BigInteger& BigInteger::operator++() {
    if (is_zero()) {
        *this = 1;
        return *this;
    }
    if (*this == -1) {
        *this = 0;
        return *this;
    }
    if (sign == plus) {
        return plus_one();
    } else {
        return minus_one();
    }
}

BigInteger BigInteger::operator++(int) {
    BigInteger answer = *this;
    ++(*this);
    return answer;
}

BigInteger& BigInteger::operator--() {
    if (this->is_zero()) {
        *this = -1;
        return *this;
    }
    if (sign == plus) {
        return minus_one();
    } else {
        return plus_one();
    }
}

BigInteger BigInteger::operator--(int) {
    BigInteger answer = *this;
    --(*this);
    return answer;
}

BigInteger::operator double() const {
    double result = 0.0;
    double base = 1.0;
    for (int i = 0; i < get_size(); ++i) {
        result += (arr[i] - '0') * base;
        base *= 10.0;
    }
    if (sign == minus) {
        result = -result;
    }
    return result;
}

BigInteger operator""_bi(const char* BigInteger_string, size_t len) {
    BigInteger answer;
    for (int i = len - 1; i >= 0; i--) {
        if (BigInteger_string[i] == '-')
            answer.change_sign();
        answer.push_back(BigInteger_string[i]);
    }
    return answer;
}

BigInteger operator""_bi(unsigned long long bigInt) {
    return BigInteger((long long)bigInt);
}
BigInteger& gcd (BigInteger& lhs, BigInteger& rhs) {
    while (lhs && rhs) {
        if (lhs <= rhs) {
            rhs %= lhs;
        } else {
            lhs %= rhs;
        }
    }
    if (!lhs)
        return rhs;
    else
        return lhs;
}

std::istream& operator>>(std::istream& in, BigInteger& bigInt) {
    bigInt.clear();
    char current = in.get();
    while ((current != '+') && (current != '-') && ((current > '9') || (current < '0'))) {
        current = in.get();
    }
    if (current == '-') {
        bigInt.change_sign();
        current = in.get();
    } else if (current == '+') {
        current = in.get();
    }
    while ((current <= '9') && (current >= '0')){
        bigInt.push_back(current);
        current = in.get();
    }
    bigInt.reverse();
    return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& bigInt) {
    out << bigInt.toString();
    return out;
}

class Rational {
private:
    BigInteger numerator, denominator;
public:
    Rational();

    Rational(int n);

    Rational(const BigInteger& n);

    Rational(const Rational& other) = default;

    Rational& operator=(const Rational& other) = default;

    ~Rational() = default;

    BigInteger get_numerator() const;

    BigInteger get_denominator() const;

    void reduce();

    Rational& operator+=(const Rational& other);

    Rational& operator-=(const Rational& other);

    Rational& operator*=(const Rational& other);

    Rational& operator/=(const Rational& other);

    bool operator==(const Rational& other) const = default;

    Rational operator-() const;

    std::string toString() const;

    std::string asDecimal(size_t) const;

    explicit operator double() const;
};

Rational::Rational() : numerator(0), denominator(1) {};

Rational::Rational(int n) : numerator(n), denominator(1) {};

Rational::Rational(const BigInteger& n) : numerator(n), denominator(1) {};

BigInteger Rational::get_numerator() const {
    return numerator;
}

BigInteger Rational::get_denominator() const {
    return denominator;
}

void Rational::reduce() {
    BigInteger numerator_abs = numerator.abs();
    BigInteger denominator_abs = denominator.abs();
    BigInteger divider = gcd (numerator_abs, denominator_abs);
    numerator /= divider;
    denominator /= divider;
    if (denominator.get_sign() == minus) {
        numerator.change_sign();
        denominator.change_sign();
    }
}

Rational& Rational::operator+=(const Rational& other) {
    numerator = numerator * other.denominator + denominator * other.numerator;
    denominator *= other.denominator;
    reduce();
    return *this;
}

Rational& Rational::operator-=(const Rational& other) {
    numerator = numerator * other.denominator - denominator * other.numerator;
    denominator *= other.denominator;
    reduce();
    return *this;
}

Rational& Rational::operator*=(const Rational& other) {
    numerator *= other.numerator;
    denominator *= other.denominator;
    reduce();
    return *this;
}

Rational& Rational::operator/=(const Rational& other) {
    numerator *= other.denominator;
    denominator *= other.numerator;
    reduce();
    return *this;
}

Rational operator+(const Rational& a, const Rational& b) {
    Rational answer = a;
    answer += b;
    return answer;
}

Rational operator-(const Rational& a, const Rational& b) {
    Rational answer = a;
    answer -= b;
    return answer;
}

Rational operator*(const Rational& a, const Rational& b) {
    Rational answer = a;
    answer *= b;
    return answer;
}

Rational operator/(const Rational& a, const Rational& b) {
    Rational answer = a;
    answer /= b;
    return answer;
}

Rational Rational::operator-() const {
    Rational answer = *this;
    if (!answer.numerator.is_zero())
        answer.numerator.change_sign();
    return answer;
}

std::strong_ordering operator<=>(const Rational& a, const Rational& b)  {
    return (a.get_numerator()* b.get_denominator()) <=> (a.get_denominator() * b.get_numerator());
}

std::string Rational::toString() const {
    int size = numerator.get_size() + 1;
    if (numerator.get_sign() == minus) {
        size++;
    }
    if (denominator != 1) {
        size += (denominator.get_size() + 1);
    }
    char answer[size];
    int now = 0;
    std::string numerator_string = numerator.toString();
    for (int i = 0; i < int(numerator_string.size()); i++) {
        answer[now++] = numerator_string[i];
    }
    if (denominator != 1) {
        std::string denominator_string = denominator.toString();
        answer[now++] = '/';
        for (int i = 0; i < denominator.get_size(); i++) {
            answer[now++] = denominator_string[i];
        }
    }
    answer[now] = '\0';
    return answer;
}

Rational::operator double() const {
    return static_cast<double>(numerator) / static_cast<double>(denominator);
}

std::string Rational::asDecimal(size_t precision = 0) const{
    if (precision == 0) {
        return (numerator / denominator).toString();
    }
    BigInteger new_numerator = numerator % denominator;
    new_numerator.multiply(precision);
    BigInteger floor = numerator / denominator;

    std::string floor_string = (floor).toString();
    std::string fraction_string = (new_numerator.abs() / denominator).toString();

    std::string answer;
    if ((floor == 0) && (numerator.get_sign() == minus)) {
        answer.push_back('-');
    }
    for (int i = 0; i < int(floor_string.size()); i++) {
        answer.push_back(floor_string[i]);
    }
    answer.push_back('.');
    for (int i = 0; i < int(precision) - int((fraction_string).size()); i++) {
        answer.push_back('0');
    }
    for (int i = 0; i < std::min(int((fraction_string).size()), int(precision)); i++) {
        answer.push_back(fraction_string[i]);
    }
    return answer;
}

std::istream& operator>>(std::istream& in, Rational& r) {
    int num;
    in >> num;
    r = num;
    return in;
}
//BigInteger finish


template<const size_t M, const size_t N, typename Field = Rational>
class Matrix {
    Field** arr;
public:
    Matrix() {
        arr = new Field*[M];
        for (size_t i = 0; i < M; i++) {
            arr[i] = new Field[N];
            std::fill(arr[i], arr[i] + N, Field(0));
        }
    }

    Matrix(const std::initializer_list<std::vector<Field>>& l) {
        arr = new Field*[M];
        size_t i = 0;
        for (const std::vector<Field>& row : l) {
            arr[i] = new Field[N];
            size_t j = 0;
            for (const Field& element : row) {
                arr[i][j++] = element;
            }
            ++i;
        }
    }

    Matrix(const Matrix& other) {
        arr = new Field*[M];
        for (size_t i = 0; i < M; i++) {
            arr[i] = new Field[N];
            std::copy(other[i], other[i] + N, arr[i]);
        }
    }

    void swap(Matrix<M, N, Field>& other) {
        std::swap(arr, other.arr);
    }

    Matrix& operator=(Matrix other) {
        swap(other);
        return *this;
    }

    Field* operator[](size_t index) {
        return arr[index];
    }

    const Field* operator[](size_t index) const {
        return arr[index];
    }

    Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other) {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                arr[i][j] += other.arr[i][j];
            }
        }
        return *this;
    }

    Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& other) const {
        Matrix<M, N, Field> answer = *this;
        answer += other;
        return answer;
    }

    Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other) {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                arr[i][j] -= other.arr[i][j];
            }
        }
        return *this;
    }

    Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& other) const {
        Matrix<M, N, Field> answer = *this;
        answer -= other;
        return answer;
    }

    Matrix<M, N, Field>& operator*=(const Field& f) {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                arr[i][j] *= f;
            }
        }
        return *this;
    }

    Matrix operator*(const Field& f) const {
        Matrix<M, N, Field> answer = *this;
        answer *= f;
        return answer;
    }

    Matrix<N, M, Field> transposed() const {
        Matrix<N, M, Field> answer;
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < M; j++) {
                answer[i][j] = arr[j][i];
            }
        }
        return answer;
    }

    template<size_t K = M>
    std::enable_if_t<K == N, Field> trace() const {
        Field answer = Field(0);
        for (size_t i = 0; i < N; i++) {
            answer += arr[i][i];
        }
        return answer;
    }

    template<size_t K = M>
    std::enable_if_t<K == N, Field> det() const {
        Field ans = Field(1);
        Matrix<M, N, Field> answer = *this;
        size_t now = 0;
        for (size_t i = 0; i < N; i++) {
            size_t num = 0;
            bool x = false;
            for (size_t j = now; j < M; j++) {
                if (answer.arr[j][i] != Field(0)) {
                    num = j;
                    x = true;
                    break;
                }
            }
            if (!x) {
                return 0;
            }
            std::swap(answer.arr[now], answer.arr[num]);
            if (num % 2 != now % 2) {
                ans *= Field(-1);
            }
            for (size_t j = now + 1; j < M; j++) {
                Field coefficient = answer.arr[j][i] / answer.arr[now][i];
                if (coefficient == Field(0))
                    continue;
                answer.arr[j][i] = Field(0);
                for (size_t k = i + 1; k < N; k++) {
                    answer.arr[j][k] -= answer.arr[now][k] * coefficient;
                }
            }
            now++;
        }
        for (size_t i = 0; i < N; i++) {
            ans *= answer[i][i];
        }
        return ans;
    }

    template<size_t K = M>
    std::enable_if_t<K == N, Matrix<M, K, Field>&> operator*=(const Matrix<N, K, Field>& other) {
        Matrix<M, K, Field> temp = *this;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < K; j++) {
                arr[i][j] = Field(0);
                for (size_t k = 0; k < N; k++) {
                    arr[i][j] += (temp[i][k] * other[k][j]);
                }
            }
        }
        return *this;
    }
    int rank() const {
        const Matrix<M, N, Field>& help = get_step_matrix();
        int answer = 0;
        for (size_t i = 0; i < M; i++) {
            bool x = false;
            for (size_t j = i; j < N; j++) {
                if (help[i][j] != Field(0)) {
                    x = true;
                    break;
                }
            }
            if (!x)
                return answer;
            answer++;
        }
        return answer;
    }

    std::vector<Field> getRow(unsigned int index) const {
        std::vector<Field> answer(N, Field(0));
        for (size_t i = 0; i < N; i++) {
            answer[i] = arr[index][i];
        }
        return answer;
    }

    std::vector<Field> getColumn(unsigned int index) const {
        std::vector<Field> answer(M, Field(0));
        for (size_t i = 0; i < M; i++) {
            answer[i] = arr[i][index];
        }
        return answer;
    }



    template<size_t K = M>
    std::enable_if_t<K == N, Matrix<M, N, Field>> inverted() const {
        Matrix<M, N, Field> answer;
        for (size_t i = 0; i < N; i++) {
            answer[i][i] = Field(1);
        }
        Matrix<M, N, Field> help = *this;
        for (size_t i = 0; i < N; i++) {
            size_t num = 0;
            for (size_t j = i; j < N; j++) {
                if (help.arr[j][i] != Field(0)) {
                    num = j;
                    break;
                }
            }
            std::swap(help.arr[i], help.arr[num]);
            std::swap(answer.arr[i], answer.arr[num]);
            for (size_t j = i + 1; j < N; j++) {
                Field coefficient = help.arr[j][i] / help.arr[i][i];
                if (coefficient == Field(0))
                    continue;
                for (size_t k = 0; k < N; k++) {
                    help.arr[j][k] -= help.arr[i][k] * coefficient;
                    answer.arr[j][k] -= answer.arr[i][k] * coefficient;
                }
            }
        }
        for (size_t i = 0; i < N; i++) {
            Field coefficient = help[i][i];
            for (size_t j = 0; j < N; j++) {
                help[i][j] /= coefficient;
                answer[i][j] /= coefficient;
            }
        }
        for (int i = N - 1; i >= 0; i--) {
            for (size_t j = i + 1; j < N; j++) {
                Field coefficient = help[i][j];
                for (size_t k = 0; k < N; k++) {
                    answer[i][k] -= answer[j][k] * coefficient;
                }
            }
        }
        return answer;
    }

    template<size_t K = M>
    std::enable_if_t<K == N, void> invert() {
        *this = this->inverted();
    }

    bool operator==(const Matrix<M, N, Field>& other) const {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                if (arr[i][j] != other[i][j])
                    return false;
            }
        }
        return true;
    }

    ~Matrix() {
        for (size_t i = 0; i < M; i++) {
            delete[] arr[i];
        }
        delete[] arr;
    }
private:
    Matrix<M, N, Field> get_step_matrix() const {
        Matrix<M, N, Field> answer = *this;
        size_t now = 0;
        for (size_t i = 0; i < N; i++) {
            size_t num = 0;
            bool x = false;
            for (size_t j = now; j < M; j++) {
                if (answer.arr[j][i] != Field(0)) {
                    num = j;
                    x = true;
                    break;
                }
            }
            if (!x) {
                continue;
            }
            std::swap(answer.arr[now], answer.arr[num]);
            for (size_t j = now + 1; j < M; j++) {
                Field coefficient = answer.arr[j][i] / answer.arr[now][i];
                if (coefficient == Field(0))
                    continue;
                answer.arr[j][i] = Field(0);
                for (size_t k = i + 1; k < N; k++) {
                    answer.arr[j][k] -= answer.arr[now][k] * coefficient;
                }
            }
            now++;
        }

        return answer;
    }
};

template<size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& a, const Matrix<N, K, Field>& b) {
    Matrix<M, K, Field> answer;
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < K; j++) {
            for (size_t k = 0; k < N; k++) {
                answer[i][j] += (a[i][k] * b[k][j]);
            }
        }
    }
    return answer;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(const Field& f, const Matrix<M, N, Field>& m) {
    Matrix<M, N, Field> answer = m;
    answer *= f;
    return answer;
}

template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field>& operator*=(const Field& f, Matrix<M, N, Field>& m) {
    m *= f;
    return m;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
