struct _Point;

struct _Point {
    double x;
    double y;

    _Point() : x(0), y(0) {}
    _Point(double x, double y) : x(x), y(y) {}
};

typedef struct _Point Point;

#define THREADS 8