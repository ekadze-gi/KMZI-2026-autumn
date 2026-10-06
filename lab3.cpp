#include <gmpxx.h>
#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

using Big = mpz_class;

struct Point {
    Big X;
    Big Y;
    Big Z;
};

static Big p;
static Big a;
static Big b;
static Big q;
static Point G;

static Big from_hex(const char* s) {
    Big x;
    x.set_str(s, 16);
    return x;
}

static std::string to_hex(const Big& x) {
    std::string s = x.get_str(16);
    for (char& c : s) {
        if (c >= 'a' && c <= 'f') c = (char)(c - 32);
    }
    if (s.size() < 64) s = std::string(64 - s.size(), '0') + s;
    return s;
}

static Big mod(const Big& x, const Big& m) {
    Big r = x % m;
    if (r < 0) r += m;
    return r;
}

static Big inv(const Big& x, const Big& m) {
    Big r;
    mpz_invert(r.get_mpz_t(), x.get_mpz_t(), m.get_mpz_t());
    return r;
}

static int bit(const Big& x, int i) {
    return mpz_tstbit(x.get_mpz_t(), i);
}

namespace constant_time {

int equal(const Big& x, const Big& y) {
    unsigned acc = 0;
    for (int i = 0; i < 256; i++) {
        int xi = bit(x, i);
        int yi = bit(y, i);
        acc |= (unsigned)(xi ^ yi);
    }
    return acc == 0;
}

int points_equal(const Point& P, const Point& Q) {
    Big z1s = mod(P.Z * P.Z, p);
    Big z2s = mod(Q.Z * Q.Z, p);
    Big u1 = mod(P.X * z2s, p);
    Big u2 = mod(Q.X * z1s, p);
    Big s1 = mod(P.Y * z2s * Q.Z, p);
    Big s2 = mod(Q.Y * z1s * P.Z, p);
    int same_x = equal(u1, u2);
    int same_y = equal(s1, s2);
    return same_x & same_y;
}

}

namespace algebra {

Big add_p(const Big& x, const Big& y) {
    return mod(x + y, p);
}

Big sub_p(const Big& x, const Big& y) {
    return mod(x - y, p);
}

Big mul_p(const Big& x, const Big& y) {
    return mod(x * y, p);
}

Big sqr_p(const Big& x) {
    return mul_p(x, x);
}

Big inv_p(const Big& x) {
    return inv(x, p);
}

Point inf() {
    Point O;
    O.X = 1;
    O.Y = 1;
    O.Z = 0;
    return O;
}

Point make_point(const Big& x, const Big& y) {
    Point P;
    P.X = x;
    P.Y = y;
    P.Z = 1;
    return P;
}

int to_affine(const Point& P, Big& x, Big& y) {
    if (P.Z == 0) return 0;
    Big zi = inv_p(P.Z);
    Big z2 = sqr_p(zi);
    Big z3 = mul_p(z2, zi);
    x = mul_p(P.X, z2);
    y = mul_p(P.Y, z3);
    return 1;
}

int on_curve(const Big& x, const Big& y) {
    Big left = sqr_p(y);
    Big right = add_p(add_p(mul_p(sqr_p(x), x), mul_p(a, x)), b);
    return left == right;
}

Point dbl(const Point& P) {
    if (P.Z == 0 || P.Y == 0) return inf();

    Big YY = sqr_p(P.Y);
    Big S = mul_p(4, mul_p(P.X, YY));
    Big Z4 = sqr_p(sqr_p(P.Z));
    Big M = add_p(mul_p(3, sqr_p(P.X)), mul_p(a, Z4));

    Point R;
    R.X = sub_p(sqr_p(M), mul_p(2, S));
    R.Y = sub_p(mul_p(M, sub_p(S, R.X)), mul_p(8, sqr_p(YY)));
    R.Z = mul_p(2, mul_p(P.Y, P.Z));
    return R;
}

Point add(const Point& P, const Point& Q) {
    if (P.Z == 0) return Q;
    if (Q.Z == 0) return P;
    if (constant_time::points_equal(P, Q)) return dbl(P);

    Big Z1Z1 = sqr_p(P.Z);
    Big Z2Z2 = sqr_p(Q.Z);
    Big U1 = mul_p(P.X, Z2Z2);
    Big U2 = mul_p(Q.X, Z1Z1);
    Big S1 = mul_p(P.Y, mul_p(Z2Z2, Q.Z));
    Big S2 = mul_p(Q.Y, mul_p(Z1Z1, P.Z));
    Big H = sub_p(U2, U1);
    Big R = sub_p(S2, S1);
    if (H == 0) return inf();

    Big H2 = sqr_p(H);
    Big H3 = mul_p(H, H2);
    Big V = mul_p(U1, H2);

    Point T;
    T.X = sub_p(sub_p(sqr_p(R), H3), mul_p(2, V));
    T.Y = sub_p(mul_p(R, sub_p(V, T.X)), mul_p(S1, H3));
    T.Z = mul_p(mul_p(H, P.Z), Q.Z);
    return T;
}

Point mul(const Big& k, const Point& P) {
    Point R = inf();
    for (int i = 255; i >= 0; i--) {
        R = dbl(R);
        if (bit(k, i)) R = add(R, P);
    }
    return R;
}

void init() {
    p = from_hex("8000000000000000000000000000000000000000000000000000000000000431");
    a = from_hex("7");
    b = from_hex("5FBFF498AA938CE739B8E022FBAFEF40563F6E6A3472FC2A514C0CE9DAE23B7E");
    q = from_hex("8000000000000000000000000000000150FE8A1892976154C59CFC193ACCF5B3");
    G = make_point(
        from_hex("2"),
        from_hex("08E2A8A0E65147D4BD6316030E16D19C85C97F0A9CA267122B96ABBCEA7E8FC8"));
}

}

namespace protocols {

Point public_key(const Big& d) {
    return algebra::mul(d, G);
}

Point shared_secret(const Big& d, const Point& Qpeer) {
    return algebra::mul(d, Qpeer);
}

int gost_sign(const Big& d, const Big& e_in, const Big& k, Big& r, Big& s) {
    Big e = mod(e_in, q);
    if (e == 0) e = 1;

    Point C = algebra::mul(k, G);
    Big x;
    Big y;
    if (!algebra::to_affine(C, x, y)) return 0;

    r = mod(x, q);
    if (r == 0) return 0;

    s = mod(r * d + k * e, q);
    if (s == 0) return 0;
    return 1;
}

int gost_verify(const Big& r, const Big& s, const Point& Qpub, const Big& e_in) {
    if (r <= 0 || s <= 0 || r >= q || s >= q) return 0;

    Big e = mod(e_in, q);
    if (e == 0) e = 1;

    Big v = inv(e, q);
    Big z1 = mod(s * v, q);
    Big z2 = mod(-(r * v), q);

    Point C = algebra::add(algebra::mul(z1, G), algebra::mul(z2, Qpub));
    Big x;
    Big y;
    if (!algebra::to_affine(C, x, y)) return 0;
    return mod(x, q) == r;
}

}

static void print_big(const char* name, const Big& x) {
    std::cout << name << " = " << to_hex(x) << "\n";
}

static void print_point(const char* name, const Point& P) {
    Big x;
    Big y;
    if (!algebra::to_affine(P, x, y)) {
        std::cout << name << " = O\n";
        return;
    }
    std::cout << name << ".x = " << to_hex(x) << "\n";
    std::cout << name << ".y = " << to_hex(y) << "\n";
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    algebra::init();

    std::cout << "Лабораторная работа 3, вариант 6\n";
    std::cout << "ГОСТ Р 34.10-2012, кривая Вейерштрасса, проективные координаты Якоби\n";
    std::cout << "Constant-time: сравнение точек XOR-аккумулятором без early exit\n\n";

    print_big("p", p);
    print_big("a", a);
    print_big("b", b);
    print_big("q", q);
    print_point("G", G);
    std::cout << "G на кривой: " << (algebra::on_curve(G.X, G.Y) ? "да" : "нет") << "\n\n";

    Big d = from_hex("7A929ADE789BB9BE10ED359DD39A72C11B60961F49397EEE1D19CE9891EC3B28");
    Point Q = protocols::public_key(d);
    Point Qref = algebra::make_point(
        from_hex("7F2B49E270DB6D90D8595BEC458B50C58585BA1D4E9B788F6689DBD8E56FD80B"),
        from_hex("26F1B489D6701DD185C8413A977B3CBBAF64D1C593D26627DFFB101A87FF77DA"));
    print_big("d", d);
    print_point("Q", Q);
    std::cout << "Q = d*G (RFC 7091): " << (constant_time::points_equal(Q, Qref) ? "да" : "нет") << "\n\n";

    Big dA = from_hex("1A2B3C4D5E6F708192A3B4C5D6E7F809112233445566778899AABBCCDDEE0011");
    Big dB = from_hex("4F00112233445566778899AABBCCDDEEFF00112233445566778899AABBCC0001");
    Point QA = protocols::public_key(dA);
    Point QB = protocols::public_key(dB);
    Point KA = protocols::shared_secret(dA, QB);
    Point KB = protocols::shared_secret(dB, QA);
    std::cout << "ECDH\n";
    print_big("dA", dA);
    print_point("QA", QA);
    print_big("dB", dB);
    print_point("QB", QB);
    print_point("KA", KA);
    print_point("KB", KB);
    int ecdh_ok = constant_time::points_equal(KA, KB);
    std::cout << "CT-сравнение KA == KB: " << (ecdh_ok ? "да" : "нет") << "\n\n";

    Big e = from_hex("2DFBC1B372D89A1188C09C52E0EEC61FCE52032AB1022E8E67ECE6672B043EE5");
    Big k = from_hex("77105C9B20BCD3122823C8CF6FCC7B956DE33814E95B7FE64FED924594DCEAB3");
    Big r;
    Big s;
    std::cout << "ГОСТ Р 34.10-2012, тестовый вектор RFC 7091\n";
    print_big("e = H(M)", e);
    print_big("k", k);
    int signed_ok = protocols::gost_sign(d, e, k, r, s);
    print_big("r", r);
    print_big("s", s);
    Big r_ref = from_hex("41AA28D2F1AB148280CD9ED56FEDA41974053554A42767B83AD043FD39DC0493");
    Big s_ref = from_hex("01456C64BA4642A1653C235A98A60249BCD6D3F746B631DF928014F6C5BF9C40");
    std::cout << "подпись построена: " << (signed_ok ? "да" : "нет") << "\n";
    std::cout << "r совпал: " << (r == r_ref ? "да" : "нет") << "\n";
    std::cout << "s совпал: " << (s == s_ref ? "да" : "нет") << "\n";
    int ver_ok = protocols::gost_verify(r, s, Q, e);
    std::cout << "проверка подписи: " << (ver_ok ? "принята" : "отвергнута") << "\n";
    Big s_bad = s ^ 1;
    print_big("s xor 1", s_bad);
    int ver_bad = protocols::gost_verify(r, s_bad, Q, e);
    std::cout << "проверка после порчи s: " << (ver_bad ? "принята" : "отвергнута") << "\n";

    int all_ok = algebra::on_curve(G.X, G.Y) && constant_time::points_equal(Q, Qref) && ecdh_ok && signed_ok &&
                 (r == r_ref) && (s == s_ref) && ver_ok && !ver_bad;
    std::cout << "\nИтог: " << (all_ok ? "PASS" : "FAIL") << "\n";
    return all_ok ? 0 : 1;
}
