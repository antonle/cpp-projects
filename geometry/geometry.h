#include <iostream>
#include <vector>
#include <cmath>

static const double EPS = 1e-9;

double abs(double number) {
    return number < 0 ? -number : number;
}

bool isEqual(double lhs, double rhs) {
    return abs(lhs - rhs) <= EPS;
}

class Point {
public:
    Point();
    Point(double x, double y, double z);
    Point(const Point& other);
    Point& operator=(const Point& other);

    double getXCoor() const;
    double getYCoor() const;
    double getZCoor() const;

private:
    double x_coor_;
    double y_coor_;
    double z_coor_;
};

Point::Point(): x_coor_(0), y_coor_(0), z_coor_(0) {}

Point::Point(double x, double y, double z): x_coor_(x), y_coor_(y), z_coor_(z) {}

Point::Point(const Point& other): x_coor_(other.x_coor_), y_coor_(other.y_coor_), z_coor_(other.z_coor_) {}

Point& Point::operator=(const Point& other) {
    x_coor_ = other.x_coor_;
    y_coor_ = other.y_coor_;
    z_coor_ = other.z_coor_;

    return *this;
}

double Point::getXCoor() const {
    return x_coor_;
}

double Point::getYCoor() const {
    return y_coor_;
}

double Point::getZCoor() const {
    return z_coor_;
}

bool operator==(const Point& lhs, const Point& rhs) {
    return isEqual(lhs.getXCoor(), rhs.getXCoor()) && isEqual(lhs.getYCoor(), rhs.getYCoor()) && isEqual(lhs.getZCoor(), rhs.getZCoor());
}

bool operator!=(const Point& lhs, const Point& rhs) {
    return !(lhs == rhs);
}

double distance(const Point& point1, const Point& point2) {
    double x1 = point1.getXCoor(), y1 = point1.getYCoor();
    double x2 = point2.getXCoor(), y2 = point2.getYCoor();

    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

Point getCenter(const Point& first, const Point& second) {
    double x1 = first.getXCoor(), y1 = first.getXCoor();
    double x2 = second.getXCoor(), y2 = second.getYCoor();

    return Point((x1 + x2) / 2, (y1 + y2) / 2, 0);
}

class Vector {
public:
    Vector();
    Vector(const Point& end);
    Vector(const Point& start, const Point& end);
    Vector(const Vector& other);
    Vector& operator=(const Vector& other);

    Point getEnd() const;

    double length() const;

    Vector& operator*=(const Vector& rhs);
    Vector& operator+=(const Vector& rhs);

private:
    Point end_;
};

Vector::Vector(): end_(0, 0, 0) {}

Vector::Vector(const Point& end): end_(end) {}

Vector::Vector(const Point& start, const Point& end): end_(end.getXCoor() - start.getXCoor(), end.getYCoor() - start.getYCoor(), end.getZCoor() - start.getZCoor()) {}

Vector::Vector(const Vector& other): end_(other.end_) {}

Vector& Vector::operator=(const Vector& other) {
    end_ = other.end_;

    return *this;
}

Point Vector::getEnd() const {
    return end_;
}

double Vector::length() const {
    return sqrt(end_.getXCoor() * end_.getXCoor() + end_.getYCoor() * end_.getYCoor() + end_.getZCoor()* end_.getZCoor());
}

Vector& Vector::operator*=(const Vector& rhs) {
    double x1 = getEnd().getXCoor(), y1 = getEnd().getYCoor(), z1 = getEnd().getZCoor();
    double x2 = rhs.getEnd().getXCoor(), y2 = rhs.getEnd().getYCoor(), z2 = rhs.getEnd().getZCoor();

    double new_x = y1 * z2 - y2 * z1;
    double new_y = x2 * z1 - x1 * z2;
    double new_z = x1 * y2 - x2 * y1;

    Point new_point(new_x, new_y, new_z);
    end_ = new_point;

    return *this;
}

Vector& Vector::operator+=(const Vector& rhs) {
    Point new_point(end_.getXCoor() + rhs.getEnd().getXCoor(), end_.getYCoor() + rhs.getEnd().getYCoor(), end_.getZCoor() + rhs.getEnd().getZCoor());
    end_ = new_point;

    return *this;
}

Vector operator*(const Vector& lhs, const Vector& rhs) {
    Vector copy = lhs;
    copy *= rhs;

    return copy;
}

Vector operator+(const Vector& lhs, const Vector& rhs) {
    Vector copy = lhs;
    copy += rhs;
print(a)

    return copy;
}

bool operator==(const Vector& lhs, const Vector& rhs) {
    return lhs.getEnd() == rhs.getEnd();
}

bool operator!=(const Vector& lhs, const Vector& rhs) {
    return !(lhs == rhs);
}

class Line {
public:
    Line();
    Line(const Point& first_point, const Point& second_point);
    Line(double koef, double shift);
    Line(const Point& point, double koef);
    Line(const Line& other);
    Line& operator=(const Line& other);

    Point getFirstPoint() const;
    Point getSecondPoint() const;

private:
    Point first_point_;
    Point second_point_;
};

Line::Line(): first_point_(0, 0, 0), second_point_(0, 0, 0) {}

Line::Line(const Point& first_point, const Point& second_point): first_point_(first_point), second_point_(second_point) {}

Line::Line(double koef, double shift): first_point_(shift, 0, 0), second_point_(1, koef + shift, 0) {}

Line::Line(const Point& point, double koef): first_point_(point), second_point_(0, point.getYCoor() - koef * point.getXCoor(), 0) {}

Line::Line(const Line& other): first_point_(other.first_point_), second_point_(other.second_point_) {}

Line& Line::operator=(const Line& other) {
    first_point_ = other.first_point_;
    second_point_ = other.second_point_;

    return *this;
}

Point Line::getFirstPoint() const {
    return first_point_;
}

Point Line::getSecondPoint() const {
    return second_point_;
}

bool operator==(const Line& lhs, const Line& rhs) {
    Vector v1(lhs.getFirstPoint(), lhs.getSecondPoint());
    Vector v2(rhs.getFirstPoint(), rhs.getSecondPoint());
    Vector v3(lhs.getFirstPoint(), rhs.getFirstPoint());

    return v1 * v2 == Vector({0, 0, 0}) && v1 * v3 == Vector({0, 0, 0});
}

bool operator!=(const Line& lhs, const Line& rhs) {
    return !(lhs == rhs);
}

class Shape {
public:
    virtual double perimeter() const;
    virtual double area() const;

    virtual bool operator==(const Shape& another) const;
    virtual bool operator!=(const Shape& another) const;
    virtual bool isCongruentTo(const Shape& another) const;
};

bool Shape::operator!=(const Shape& another) const {
    return !(*this == another);
}

class Polygon: public Shape {
public:
    Polygon() = default;
    Polygon(std::vector<Point>& points);
    Polygon(std::initializer_list<Point>& points);

    template <typename T = Point, typename... Args>
    Polygon(T point, Args... points);

    size_t verticesCount() const;
    std::vector<Point> getVertices() const;

    bool isConvex() const;

    double perimeter() const override;
    double area() const override;

    bool operator==(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;

protected:
    std::vector<Point> points_;
};

Polygon::Polygon(std::vector<Point>& points): points_(points) {}

Polygon::Polygon(std::initializer_list<Point>& points): points_(points) {}

template <typename T, typename... Args>
Polygon::Polygon(T point, Args... points): Polygon(points...) {
    points_.push_back(point);
}

size_t Polygon::verticesCount() const {
    return points_.size();
}

std::vector<Point> Polygon::getVertices() const {
    return points_;
}

bool Polygon::isConvex() const {
    bool flag;
    double copy;
    for (size_t i = 0; i < points_.size() - 2; ++i) {
        Vector first(points_[i], points_[i + 1]);
        Vector second(points_[i + 1], points_[i + 2]);
        Vector proc = first * second;

        flag = i == 0 || copy * proc.getEnd().getZCoor();
        copy = proc.getEnd().getZCoor();
    }

    return flag;
}

double Polygon::perimeter() const {
    double ans = 0;

    for (size_t i = 0; i < points_.size() - 1; ++i) {
        Point first = points_[i], second = points_[i + 1];
        ans += distance(first, second);
    }

    return ans;
}

double Polygon::area() const {
    Vector ans = Vector(points_[0], points_[1]) * Vector(points_[0], points_[2]);

    for (size_t i = 2; i < points_.size(); ++i) {
        ans += Vector(points_[0], points_[i]) * Vector(points_[0], points_[i + 1]);
    }

    return ans.length();
}

bool Polygon::operator==(const Shape& another) const {
    const Polygon* other = static_cast<const Polygon*>(&another);

    bool flag = false;
    for (size_t i = 0; i < points_.size(); ++i) {
        for (size_t j = 0; j < other->points_.size(); ++j) {
            if (points_[i] == other->points_[i]) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            return false;
        }
        flag = false;
    }

    return true;
}

bool Polygon::isCongruentTo(const Shape& another) const {

}

class Eclipse: public Shape {
public:
    Eclipse();
    Eclipse(const Point& first_focus, const Point& second_focus, double distance);

    std::pair<Point, Point> focuses() const;
    Point center() const;
    double eccentricity() const;
    std::pair<Line, Line> directrices() const;

    double perimeter() const override;
    double area() const override;

    bool operator==(const Shape& another) const override;
    bool isCongruentTo(const Shape& another) const override;

protected:
    std::pair<Point, Point> focuses_;
    double distance_of_focuses_;
};

Eclipse::Eclipse(): focuses_({{0, 0, 0}, {0, 0 ,0}}), distance_of_focuses_(0) {}

Eclipse::Eclipse(const Point& first_focus, const Point& second_focus, double distance): focuses_({first_focus, second_focus}),
                                                                                         distance_of_focuses_(distance)
{}

std::pair<Point, Point> Eclipse::focuses() const {
    return focuses_;
}

Point Eclipse::center() const {
    return Point((focuses_.first.getXCoor() + focuses_.second.getXCoor()) / 2,
                 (focuses_.first.getYCoor() + focuses_.second.getYCoor()) / 2,
                 0);
}

double Eclipse::eccentricity() const {
    Point cent = center();
    return 2 * (focuses_.second.getXCoor() - cent.getXCoor()) / distance_of_focuses_;
}

std::pair<Line, Line> Eclipse::directrices() const {
    Point cent = center();
    double x_coor = cent.getXCoor() + distance_of_focuses_ / 2 * eccentricity();
    Line first_directrice({x_coor, 0, 0}, {x_coor, 1, 0});
    Line second_directrice({-x_coor, 0, 0}, {-x_coor, 1, 0});

    return {first_directrice, second_directrice};
}

double Eclipse::perimeter() const {
    double e = eccentricity();
    double a = distance_of_focuses_ / 2;
    double b = a * sqrt(1 - e * e);
    double h = (a - b) * (a - b) / ((a + b) * (a + b));

    return M_PI * (a + b) * (1 + 3 * h / (10 + sqrt(4 - 3 * h)));
}

double Eclipse::area() const {
    double a = distance_of_focuses_ / 2;
    double e = eccentricity();
    double b = a * sqrt(1 - e * e);

    return M_PI * a * b;
}

bool Eclipse::operator==(const Shape& another) const {
    const Eclipse* other = static_cast<const Eclipse*>(&another);

    return focuses_ == other->focuses_ && isEqual(distance_of_focuses_, other->distance_of_focuses_) && directrices() == other->directrices();
}

bool Eclipse::isCongruentTo(const Shape& another) const {
    const Eclipse* other = static_cast<const Eclipse*>(&another);

    return focuses_ == other->focuses_ && isEqual(distance_of_focuses_, other->distance_of_focuses_);
}

class Circle: public Eclipse {
public:
    Circle();
    Circle(const Point& center, double radius);

    double radius() const;

    double perimeter() const override;
    double area() const override;
};

Circle::Circle(): Eclipse() {}

Circle::Circle(const Point& center, double radius): Eclipse(center, center, 2 * radius) {}

double Circle::radius() const {
    return distance_of_focuses_ / 2;
}

double Circle::perimeter() const {
    return 2 * M_PI * radius();
}

double Circle::area() const {
    return M_PI * radius() * radius();
}

class Rectangle: public Polygon {
public:
    Rectangle(const Point& first_point, const Point& second_point, double ratio);

    Point center() const;
    std::pair<Line, Line> diagonals() const;

    double perimeter() const override;
    double area() const override;
};

Rectangle::Rectangle(const Point& first_point, const Point& second_point, double ratio) {
    double x1 = first_point.getXCoor(), y1 = second_point.getYCoor();
    double x3 = first_point.getXCoor(), y3 = second_point.getYCoor();
    double x2 = (x1 + ratio * ratio * x3 - ratio * (y3 - y1)) / (ratio * ratio + 1);
    double y2 = (ratio * y3 - x1 + x2) / ratio;
    double x4 = x3 - x2 + x1;
    double y4 = y3 - y2 + y1;

    points_.push_back(first_point);
    points_.push_back(Point(x2, y2, 0));
    points_.push_back(second_point);
    points_.push_back(Point(x4, y4, 0));
}

Point Rectangle::center() const {
    return Point((points_[2].getXCoor() - points_[0].getXCoor()) / 2, (points_[2].getYCoor() - points_[0].getYCoor()) / 2, 0);
}

std::pair<Line, Line> Rectangle::diagonals() const {
    Line first_diag(points_[0], points_[2]);
    Line second_diag(points_[1], points_[3]);

    return {first_diag, second_diag};
}

double Rectangle::perimeter() const {
    double ab = distance(points_[0], points_[1]);
    double bc = distance(points_[1], points_[2]);

    return 2 * (ab + bc);
}

double Rectangle::area() const {
    double ab = distance(points_[0], points_[1]);
    double bc = distance(points_[1], points_[2]);

    return ab * bc;
}

class Square: public Rectangle {
public:
    Square(const Point& first_point, const Point& second_point);

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Square::Square(const Point& first_point, const Point& second_point): Rectangle(first_point, second_point, 1) {}

Circle Square::circumscribedCircle() const {
    return Circle(center(), distance(points_[0], points_[2]) / 2);
}

Circle Square::inscribedCircle() const {
    return Circle(center(), distance(points_[0], points_[1]) / 2);
}

class Triangle: public Polygon {
public:
    Triangle(const Point& first_point, const Point& second_point, const Point& third_point);

    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;

    double perimeter() const override;
    double area() const override;
};

Triangle::Triangle(const Point& first_point, const Point& second_point, const Point& third_point): Polygon(first_point, second_point, third_point) {}

Circle Triangle::circumscribedCircle() const {
    Point a = points_[0], b = points_[1], c = points_[2];

    Point ma = getCenter(b, c), mb = getCenter(a, c), mc = getCenter(a, b);

    Triangle new_tri(ma, mb, mc);
    Point o = new_tri.orthocenter();

    return Circle(o, distance(o, a));
}

Circle Triangle::inscribedCircle() const {
    Point a = points_[0], b = points_[1], c = points_[2];

    double xa = a.getXCoor(), ya = a.getYCoor();
    double xb = b.getXCoor(), yb = b.getYCoor();
    double xc = c.getXCoor(), yc = c.getYCoor();

    double ab = distance(a, b);
    double bc = distance(b, c);
    double ca = distance(c, a);

    double k = ab / ca;

    double xl = (k * xc + xb) / (1 + k), yl = (k * yc + yb) / (1 + k);
    k = ab / (k * ca);

    double xi = (k * xl + xa) / (1 + k), yi = (k * yl + ya) / (1 + k);
    Point i(xi, yi, 0);

    Vector x(a, b);
    Vector y(a, i);
    double h = ((x * y).length()) / (4 * x.length());

    return Circle(i, h);
}

Point Triangle::centroid() const {
    Point a = points_[0];
    Point b = points_[1];
    Point c = points_[2];

    Point ma = getCenter(b, c);
    return Point((2 * ma.getXCoor() + a.getXCoor()) / 3, (2 * ma.getYCoor() + a.getYCoor()) / 3, 0);
}

Point Triangle::orthocenter() const {
    Point a = points_[0];
    Point b = points_[1];
    Point c = points_[2];

    Circle out = circumscribedCircle();
    Point o = out.center();

    double xa = a.getXCoor(), ya = a.getYCoor();
    double xb = b.getXCoor(), yb = b.getYCoor();
    double xc = c.getXCoor(), yc = c.getYCoor();
    double xo = o.getXCoor(), yo = o.getYCoor();

    return Point(xa + xb + xc - 2 * xo, ya + yb + yc - 2 * yo, 0);
}

Line Triangle::EulerLine() const {
    return Line(centroid(), orthocenter());
}

Circle Triangle::ninePointsCircle() const {
    Point a = points_[0];
    Point b = points_[1];
    Point c = points_[2];

    Point ma = getCenter(b, c);
    Point mb = getCenter(a, c);
    Point mc = getCenter(a, b);

    Triangle new_tri(ma, mb, mc);

    return new_tri.circumscribedCircle();
}

double Triangle::perimeter() const {
    return distance(points_[0], points_[1]) + distance(points_[1], points_[2]), distance(points_[2], points_[0]);
}

double Triangle::area() const {
    return perimeter() * inscribedCircle().radius() / 2;
}
