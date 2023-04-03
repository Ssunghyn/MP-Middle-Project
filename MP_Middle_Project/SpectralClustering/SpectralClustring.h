struct _Point;

struct _Point {
    double x;
    double y;

    _Point(double x, double y) : x(x), y(y) {}
};

typedef struct _Point Point;