#ifndef _GRIDFNC_H_
#define _GRIDFNC_H_

#include <iostream>
#include <valarray>
#include <complex>        // should come before fftw3.h so that fftw can use it
#include "fftw3.h"
//
//    為了讓 nbplot 畫 gridfnc 必需知道 private 的 x, y
//  所以想把 nbplot 定成 friend
//  互相 include 會造成 recursive problem 而編譯不過
//  如果一定要用，就要用 forward declaration
//  不過 friend template 還有別的問題，暫時放棄這個作法
//
//  #include "tablemaker.h"  
//  #include "col.h"
//
//  class gridfnc;
//  class tablemaker
//  {
//  public:
//      void nbprint( const gridfnc<>& z );
//  };
//

using namespace std;

class gridfnc
{
protected: 
    valarray<double> x, y;
    double a, b, h;
    int n;                                                            // number of grid points = n + 1
    static int order;
    static bool trap_conv;                                            // default value = 0; use midpt rule

    // for convolution
    valarray<double> repeatL (const valarray<double>& a);        // 在左邊抄 n-1 次的 a[0], n = a.size()
    valarray<double> repeatR (const valarray<double>& a);        // 在右邊抄 n-1 次的 a[n-1]，這兩個梯形法 convol 積分要用
    valarray<double> padR (const valarray<double>& a);        // 在右邊抄 n-1 次 0, for linear convolution
    valarray<double> truncate (const valarray<double>& a);    // 砍左右得到 convol 積分結果
    
public:
    gridfnc (){}
    gridfnc (const double& a, const double& b, const int& n, double (*f)(double));
    gridfnc (const double& a, const double& b, const int& n, const double& val);    // indicator 1_{x > val}
    gridfnc (const double& a, const double& b, const int& n);                    // f(x) = x 
    gridfnc (const gridfnc& f);                                        // copy constructor
    ~gridfnc (void){}
    
    static void set_order(int ord);
    static void use_trap_conv(bool flag=1);                            // use_trap_conv(): use trap rule for convolution
                                                                    // use_trap_conv(0): use midpt rule
    
    /*    ////////////////////////////////// 目前只有 2 階版本的成員函數 //////////////////////////////////
    double integral ();                                                // Trapezoidal rule
    double inner_prod (const gridfnc& g);                            // integrate f*g from a to b using Trapezoidal rule, 假設 f.x = g.x
    double inner_prod (double (*g)(double));
    
    gridfnc derivative ();                                    // 1st derivative (function), approximated by 2nd order finite difference
    double antiderivative ( const double& x );
    gridfnc antiderivative ();

    void convolve_with (const gridfnc& g);
    void convolve_with (double (*g)(double));
    
    ////////////////////////////////// 已經有 4 階版本的 //////////////////////////////////
    double operator() ( const double& val )
    double derivative ( const double& val )
    double derivative_2nd ( const double& val )
    */
    
    // set_gridpt 改變 x，a, b, n, h, y 也跟著變
    void set_gridpt (const double& a, const double& b, const int& n);            // set x grid pt, also used in constructors
    void set_gridpt (const gridfnc& g);                            // set *this.x = g.x
    double integral ();                                                    // Trapezoidal rule
    double sup ();
    double inf ();
    double inner_prod (const gridfnc& g);                                // integrate f*g from a to b using Trapezoidal rule, 假設 f.x = g.x
    double inner_prod (double (*g)(double));
    gridfnc apply (double (*g)(double));                                    // g(*this) 反序合成，不需要 input gridfnc 的版本因為呼叫 () 就好了
    
    gridfnc derivative ();                                        // 1st derivative (function), approximated by 2nd order finite difference
    double derivative ( const double& x );                                    // 1st derivative at x approximated by interpolating quadratic function
    double derivative_2nd ( const double& x );                                // 2nd derivative at x approximated by interpolating cubic function
    // gridfnc antiderivative ();
    // double antiderivative ( const double& x );

    double root_near ( const double& x0 );                                    // Newton's method
    double root_between (double a, double b);                                        // bisection 
    double geta ();
    double getb ();
    double geth ();
    int getn ();
    int getindex(const double& val);                                        // find index i such that val is in (x[i], x[i+1])
    double getxi (const double& val);                                            // 回傳 val 左邊最近的格點
    double getyi (const double& val);                                            // 回傳 val 左邊最近的格點上的函數值
    valarray<double> getx ();
    valarray<double> gety ();
    void sety ( const valarray<double>& fncval );
    void setyi (int i, const double& val);

    // 假設兩函數用一樣的 x grid, 假設 a<0, b>0
    // 目前只能跑 double = double，如果要改 long double，要呼叫別的 fftw functions，如果要改 float，至少丟進去的 array 要先轉型成 double
    void convolve_with (const gridfnc& g);
    void convolve_with (double (*g)(double));

    // operators
    gridfnc operator+ () const;
    gridfnc operator- () const;

    double operator() (const double& val);
    gridfnc operator() (gridfnc& f);            // 超載 () 時用 const 編譯不過。可以自己合成自己。g(h) 時用的是 g 的 x grid
    gridfnc operator() (double (*f)(double));

    gridfnc& operator= (const gridfnc& f);
    gridfnc& operator= (const double& val);
    gridfnc& operator= (double (*f)(double));                // 不能丟內建函數，要自己另定義 double f(double) 才能用

    gridfnc& operator*= (const gridfnc& rhs);
    gridfnc& operator/= (const gridfnc& rhs);
    gridfnc& operator+= (const gridfnc& rhs);
    gridfnc& operator-= (const gridfnc& rhs);

    gridfnc& operator*= (const double& val);
    gridfnc& operator/= (const double& val);
    gridfnc& operator+= (const double& val);
    gridfnc& operator-= (const double& val);

    gridfnc& operator*= (double (*f)(double));
    gridfnc& operator/= (double (*f)(double));
    gridfnc& operator+= (double (*f)(double));
    gridfnc& operator-= (double (*f)(double));

    // nonmember functions                            operator+-*/ (double (*f)(double), const gridfnc& rhs) 這種的四個都不能用
    //                                                神奇的是 max (double (*f)(double), const gridfnc& rhs) 卻可以用
    //
    friend gridfnc operator* (const gridfnc& lhs, const gridfnc& rhs);
    friend gridfnc operator* (const double& val, const gridfnc& rhs);
    friend gridfnc operator* (const gridfnc& lhs, const double& val);
    friend gridfnc operator* (double (*f)(double), const gridfnc& rhs);
    friend gridfnc operator* (const gridfnc& lhs, double (*f)(double));

    friend gridfnc operator/ (const gridfnc& lhs, const gridfnc& rhs);
    friend gridfnc operator/ (const double& val, const gridfnc& rhs);
    friend gridfnc operator/ (const gridfnc& lhs, const double& val);
    friend gridfnc operator/ (double (*f)(double), const gridfnc& rhs);
    friend gridfnc operator/ (const gridfnc& lhs, double (*f)(double));

    friend gridfnc operator+ (const gridfnc& lhs, const gridfnc& rhs);
    friend gridfnc operator+ (const double& val, const gridfnc& rhs);
    friend gridfnc operator+ (const gridfnc& lhs, const double& val);
    friend gridfnc operator+ (double (*f)(double), const gridfnc& rhs);
    friend gridfnc operator+ (const gridfnc& lhs, double (*f)(double));

    friend gridfnc operator- (const gridfnc& lhs, const gridfnc& rhs);
    friend gridfnc operator- (const double& val, const gridfnc& rhs);
    friend gridfnc operator- (const gridfnc& lhs, const double& val);
    friend gridfnc operator- (double (*f)(double), const gridfnc& rhs);
    friend gridfnc operator- (const gridfnc& lhs, double (*f)(double));

    friend gridfnc max (const gridfnc& lhs, const gridfnc& rhs);    // 假設 f.x = g.x
    friend gridfnc max (const double& val, const gridfnc& rhs);
    friend gridfnc max (const gridfnc& lhs, const double& val);
    friend gridfnc max (double (*f)(double), const gridfnc& rhs);
    friend gridfnc max (const gridfnc& lhs, double (*f)(double));

    friend gridfnc min (const gridfnc& lhs, const gridfnc& rhs);
    friend gridfnc min (const double& val, const gridfnc& rhs);
    friend gridfnc min (const gridfnc& lhs, const double& val);
    friend gridfnc min (double (*f)(double), const gridfnc& rhs);
    friend gridfnc min (const gridfnc& lhs, double (*f)(double));

    friend gridfnc convolve (const gridfnc& f, const gridfnc& g);    // 用 convolve_with，呼叫者會被改變；用 convolve 不會
    friend gridfnc convolve (const gridfnc& f, double (*g)(double));
    friend gridfnc convolve (double (*g)(double), const gridfnc& f);

    // overloaded cmath functions，缺 atan2
    friend gridfnc abs (const gridfnc& f);
    friend gridfnc acos (const gridfnc& f);
    friend gridfnc asin (const gridfnc& f);
    friend gridfnc atan (const gridfnc& f);
    friend gridfnc cos (const gridfnc& f);
    friend gridfnc cosh (const gridfnc& f);
    friend gridfnc exp (const gridfnc& f);
    friend gridfnc log (const gridfnc& f);
    friend gridfnc log10 (const gridfnc& f);
    friend gridfnc sin (const gridfnc& f);
    friend gridfnc sinh (const gridfnc& f);
    friend gridfnc sqrt (const gridfnc& f);
    friend gridfnc tan (const gridfnc& f);
    friend gridfnc tanh (const gridfnc& f);
    friend gridfnc pow (const gridfnc& f, const gridfnc& g);        // 假設 f.x = g.x
    friend gridfnc pow (const gridfnc& f, const double& exponent);
    friend gridfnc pow (const double& base, const gridfnc& f);

    // overloaded special functions
    friend gridfnc ncdf (const gridfnc& f);
    friend gridfnc npdf (const gridfnc& f);
};

template<class T> int gridfnc::order = 2;
template<class T> bool gridfnc::trap_conv = 0;

// private functions for convolution
valarray<double> gridfnc::repeatL(const valarray<double>& a)
{    
    int n = a.size();
    valarray<double> result(2*n-1);
    
    for(int i=0 ; i<n-1 ; i++)
        result[i] = a[0];
    for(int i=n-1 ; i<2*n-1 ; i++)
        result[i] = a[i-n+1];

    return result;
}
valarray<double> gridfnc::repeatR(const valarray<double>& a)
{
    int n = a.size();
    valarray<double> result(2*n-1);
    
    for(int i=0 ; i<n-1 ; i++)
        result[i] = a[i];
    for(int i=n-1 ; i<2*n-1 ; i++)
        result[i] = a[n-1];

    return result;
}
valarray<double> gridfnc::padR(const valarray<double>& a)
{
    int n = a.size();
    valarray<double> result(2*n-1);
    
    for( int i=0 ; i<a.size() ; i++ )
        result[i] = a[i];

    return result;
}
valarray<double> gridfnc::truncate(const valarray<double>& a)
{
    // a.size() is assumed to be odd
    // n is the size of the original array before padding 

    int n = (a.size() + 1)/2;
    valarray<double> result(n);
    
    if( n%2 == 1 ){            // original size of the array is odd
        for(int i=0 ; i<n ; i++)
            result[i] = a[(n-1)/2 + i];
    }else{    // even; won't get the integrals at grid points; interpolation needed
        for(int i=0 ; i<n ; i++)
            result[i] = (a[n/2 - 1 + i] + a[n/2 + i])/2;
    }
    return result;
}

// constructors
gridfnc::gridfnc(const double& a, const double& b, const int& n, double (*f)(double))
{
    set_gridpt(a, b, n);
    y = x.apply(f);
}
gridfnc::gridfnc(const double& a, const double& b, const int& n, const double& val)
{
    set_gridpt(a, b, n);
    // y = val; constant function 版本，後來改為 indicator 1_{x > val}

    for( int i=0 ; i<=n ; i++ )
        y[i] = (x[i]>val) ? 1 : 0;
}
gridfnc::gridfnc(const double& a, const double& b, const int& n)
{
    set_gridpt(a, b, n);
    y = x;
}
gridfnc::gridfnc(const gridfnc& f)
{
    // f 的所有東西都是 private，可是拿的到
    
    set_gridpt(f.a, f.b, f.n);
    y = f.y;
}

void gridfnc::set_order(int ord)
{
    order = ord;
}
void gridfnc::use_trap_conv(bool flag=1)
{
    trap_conv = flag;
}

void gridfnc::set_gridpt (const double& lb, const double& ub, const int& num)
{
    a = lb;
    b = ub;
    n = num;
    h = (b-a)/n;
    
    if( x.size()!=0 )    // 本來就有 x, y 了。會呼叫 operator() 所以 x 和 y 是不能一邊算一邊蓋掉的
    {    
        valarray<double> tmpx(n+1);
        valarray<double> tmpy(n+1);
    
        tmpx[0] = a;
        tmpy[0] = (*this)(a);

        for(int i=1 ; i<=n ; i++){
            tmpx[i] = tmpx[i-1] + h;
            tmpy[i] = (*this)(tmpx[i]);
        }
        x = tmpx;
        y = tmpy;
    }
    else                // 本來沒有 x, y。
    {
        x.resize(n+1);
        y.resize(n+1);
        
        x[0] = a;
        for(int i=1 ; i<=n ; i++)
            x[i] = x[i-1] + h;
    }
}
void gridfnc::set_gridpt (const gridfnc& g)
{
    valarray<double> tmpy(g.n + 1);

    for(int i=0 ; i<=g.n ; i++)
        tmpy[i] = (*this)(g.x[i]);
    
    x = g.x;
    y = tmpy;

    // D 了一晚才找出的 bug：所有資料都要改，包括 a, b, h, n
    a = g.a;
    b = g.b;
    n = g.n;
    h = g.h;
}                            
gridfnc::integral()
{
    return h*(y.sum() - y[0]/2 - y[n]/2);
}
gridfnc::sup()
{
    // avoid discontinuity caused by numerical err at a, b
    valarray<double> tmpy = y;
    
    tmpy[0] = tmpy[1];
    tmpy[n] = tmpy[n-1];
    
    return tmpy.max();
}
gridfnc::inf()
{
    // avoid discontinuity caused by numerical err at a, b
    valarray<double> tmpy = y;
    
    tmpy[0] = tmpy[1];
    tmpy[n] = tmpy[n-1];

    return tmpy.min();
}
gridfnc::inner_prod(const gridfnc& g)
{
    return h*(((g.y)*y).sum() - (g.y[0])*y[0]/2 - (g.y[n])*y[n]/2);
}
gridfnc::inner_prod(double (*g)(double))
{
    return h*((x.apply(g)*y).sum() - (g(a))*y[0]/2 - (g(b))*y[n]/2);
}
gridfnc gridfnc::apply(double (*g)(double))
{
    gridfnc result; 

    result.set_gridpt(a, b, n);
    result.y = y.apply(g);

    return result;
}

gridfnc::derivative( const double& val )
{
    if( val < a || val > b )    return 0;
    else
    {
        int i = getindex( val );
        
        if( order==2 )
        {
            if( i==0 )
                return (-(h*(3*y[0] - 4*y[1] + y[2])) + 2*(val - x[0])*(y[0] - 2*y[1] + y[2]))/(2*h*h);
            else
                return (h*(y[i+1] - y[i-1]) + 2*val*(-2*y[i] + y[i+1] + y[i-1]) + x[i]*(4*y[i] - 2*(y[i+1] + y[i-1])))/(2*h*h);
        }
        else if( order==4 )
        {
            if( i==0 || i==1 )            // 只準到三階，還有數值穩定度的問題，因為除以 h^4，i=n 附近還沒寫
                return (-3*(val - x[0])*(val - x[0])*(y[0] - 3*y[1] + 3*y[2] - y[3]) + 6*h*(val - x[0])*(2*y[0] - 5*y[1] + 4*y[2] - y[3]) + h*h*(-11*y[0] + 18*y[1] - 9*y[2] + 2*y[3]))/(6*h*h*h);
            else{
                // 先算四個格點上的 derivative, 再用 cubic interpolation。derivative 公式來自 http://en.wikipedia.org/wiki/Five-point_stencil

                double fpp[] = { (y[i-3] - 8*y[i-2] + 8*y[i]   - y[i+1])/(12*h), 
                            (y[i-2] - 8*y[i-1] + 8*y[i+1] - y[i+2])/(12*h), 
                            (y[i-1] - 8*y[i]   + 8*y[i+2] - y[i+3])/(12*h), 
                            (y[i]   - 8*y[i+1] + 8*y[i+3] - y[i+4])/(12*h)    };

                return BarycentricCubic( &x[i-1], fpp, val );
            }
            // 只準到三階，還有數值穩定度的問題，因為除以 h^4
            // return (-3*(val - x[i-2])*(val - x[i-2])*(y[i-2] - 3*y[i-1] + 3*y[i] - y[i+1]) + 6*h*(val - x[i-2])*(2*y[i-2] - 5*y[i-1] + 4*y[i] - y[i+1]) + h*h*(-11*y[i-2] + 18*y[i-1] - 9*y[i] + 2*y[i+1]))/(6*h*h*h);
        }
        else
        {
            cout << "order is neither 2 nor 4" << endl;
            return 0;
        }
    }
}
gridfnc::derivative_2nd( const double& val )
{
    if( val < a || val > b )    return 0;
    else
    {
        int i = getindex( val );

        if( order==2 )
        {
            if( i==0 || i==1 )
                return (-((val - x[0])*(y[0] - 3*y[1] + 3*y[2] - y[3])) + h*(2*y[0] - 5*y[1] + 4*y[2] - y[3]))/(h*h*h);
            else
                return (-((val - x[i-2])*(y[i-2] - 3*y[i-1] + 3*y[i] - y[i+1])) + h*(2*y[i-2] - 5*y[i-1] + 4*y[i] - y[i+1]))/(h*h*h);
        }
        else if( order==4 )
        {
            if( i<=2 )                    // 只準到三階，還有數值穩定度的問題，因為除以 h^4，i=n 附近還沒寫
                return ( 6*(val - x[0])*(val - x[0])*(y[0] - 4*y[1] + 6*y[2] - 4*y[3] + y[4]) 
                        - 6*h*(val - x[0])*(5*y[0] - 18*y[1] + 24*y[2] - 14*y[3] + 3*y[4]) 
                        + h*h*(35*y[0] - 104*y[1] + 114*y[2] - 56*y[3] + 11*y[4]))/(12*h*h*h*h);
            else{
                // 先算四個格點上的 2nd derivative, 再用 cubic interpolation。derivative 公式來自 http://en.wikipedia.org/wiki/Five-point_stencil

                double fpp[] = {    (-y[i-3] + 16*y[i-2] -30*y[i-1] + 16*y[i]   - y[i+1])/(12*h*h), 
                            (-y[i-2] + 16*y[i-1] -30*y[i]   + 16*y[i+1] - y[i+2])/(12*h*h), 
                            (-y[i-1] + 16*y[i]   -30*y[i+1] + 16*y[i+2] - y[i+3])/(12*h*h), 
                            (-y[i]   + 16*y[i+1] -30*y[i+2] + 16*y[i+3] - y[i+4])/(12*h*h)    };

                return BarycentricCubic( &x[i-1], fpp, val );
            }    
            /* // 只準到三階，還有數值穩定度的問題，因為除以 h^4
                return ( 6*(val - x[i-3])*(val - x[i-3])*(y[i-3] - 4*y[i-2] + 6*y[i-1] - 4*y[i] + y[i+1]) 
                        - 6*h*(val - x[i-3])*(5*y[i-3] - 18*y[i-2] + 24*y[i-1] - 14*y[i] + 3*y[i+1]) 
                        + h*h*(35*y[i-3] - 104*y[i-2] + 114*y[i-1] - 56*y[i] + 11*y[i+1]))/(12*h*h*h*h);
            */
        }
        else
        {
            cout << "order is neither 2 nor 4" << endl;
            return 0;
        }
    }
}
gridfnc gridfnc::derivative()
{
    gridfnc result;

    result.set_gridpt(a, b, n);

    for(int i=1 ; i<n ; i++)
        result.y[i] = (y[i+1] - y[i-1])/(2*h);

    result.y[0] = (-3*y[0] + 4*y[1] - y[2])/(2*h);
    result.y[n] = (3*y[n] - 4*y[n-1] + y[n-2])/(2*h);

    return result;
}

gridfnc::root_near( const double& x0 )
{
    double tmpx, x = x0;
    double absTol = 0.000000000001;    // 10^{-12}
    double relTol = 0.000000000001;    // 10^{-12}

    //gridfnc fp = derivative();

    for(int i=0 ; i<100 ; i++)
    {
        tmpx = x;            // need last x for increment test 
        // x = x - (*this)(x)/fp(x);
        double fp = derivative(x);
        x = x - (*this)(x)/fp;
        if( sign(abs(x-tmpx)-(absTol + relTol*abs(x))) < 0 )
            break;
    }

    return x;
}
gridfnc::root_between (double a, double b)
{
    double x;
    double absTol = 0.000000000001;    // 10^{-12}
    double relTol = 0.000000000001;    // 10^{-12}
    int S;

    S = floor(log((b-a)/(absTol + relTol*((abs(a)<abs(b)) ? abs(a) : abs(b))))/log(2.0));
    x = (a+b)/2.0;

    for(int k=0 ; k<S ; k++){
        if( sign((*this)(x))*sign((*this)(a)) < 0 ){
            b = x;
        }else{
            a = x;
        }
        x = a + (b-a)/2.0;
    }
    return x;
}

int gridfnc::getindex(const double& val)
{
    return (int)((val - x[0])/h);
}
gridfnc::geta()
{
    return a;
}
gridfnc::getb()
{
    return b;
}
gridfnc::geth()
{
    return h;
}
gridfnc::getxi(const double& val)
{
    return x[ getindex(val) ];
}
gridfnc::getyi(const double& val)
{
    return y[ getindex(val) ];
}
int gridfnc::getn()
{
    return n;
}
valarray<double> gridfnc::getx ()
{
    return x;
}
valarray<double> gridfnc::gety ()
{
    return y;
}
void gridfnc::sety ( const valarray<double>& fncval )
{
    y = fncval;
}
void gridfnc::setyi ( int i, const double& val )
{
    y[i] = val;
}

// convolution
void gridfnc::convolve_with(const gridfnc& g)
{
    // 目前只能用 double
    // 假設 *this 和 g 的 x grid 相同
    // 假設 a<0, b>0
    // set up valarrays for test

    valarray<double> a = padR(y);                    // 在原 sequence 後面填 0
    valarray<double> b = padR(g.y);
    valarray< complex<double> > fa(2*n-1);            // 用來存 Fourier[a]
    valarray< complex<double> > fb(2*n-1);
    
    double *pa = &a[0];                                // pa -> a[0]
    double *pb = &b[0];
    fftw_complex *pfa = (fftw_complex*)&fa[0];        // pfa -> fa[0]
    fftw_complex *pfb = (fftw_complex*)&fb[0];

    // invF( F(a) * F(b) )
    fftw_execute ( fftw_plan_dft_r2c_1d ( 2*n-1, pa, pfa, FFTW_ESTIMATE ) );
    fftw_execute ( fftw_plan_dft_r2c_1d ( 2*n-1, pb, pfb, FFTW_ESTIMATE ) );
    fa *= fb;
    fftw_execute ( fftw_plan_dft_c2r_1d ( 2*n-1, pfa, pa, FFTW_ESTIMATE ) );
    
    a *= (h/(2*n-1));

    // 梯形法左右校正項，沒有這行就是中點矩形法從 a-h/2 積到 b+h/2
    if( trap_conv )    
        a -= 0.5*(repeatR(y)*repeatL(g.y) + repeatL(y)*repeatR(g.y))*h;

    // truncation & interpolation for even-length array
    y = truncate(a);
}
void gridfnc::convolve_with(double (*g)(double))
{
    gridfnc tmp(a, b, n, g);
    (*this).convolve_with(tmp);
}

// operators
gridfnc gridfnc::operator+() const
{
    return *this;
}
gridfnc gridfnc::operator-() const
{
    gridfnc result;

    result.set_gridpt(a, b, n);
    result.y = -y;
    return result;
}

gridfnc::operator()(const double& val)
{
    if( val < a || b < val){
        return 0;
    }else{
        int i = getindex(val);

        if( order == 2 )
            return y[i] + (val - x[i])*(y[i+1] - y[i])/(x[i+1] - x[i]);
        else if( order == 4 ){
            if        ( i==0 )    return BarycentricCubic( &x[0],   &y[0],   val );
            else if    ( i>=n-1 )    return BarycentricCubic( &x[n-3], &y[n-3], val );
            else                return BarycentricCubic( &x[i-1], &y[i-1], val );
        }else{
            cout << "order is neither 2 nor 4" << endl;
            return 0;
        }
    }
}
gridfnc gridfnc::operator()(gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(a, b, n);

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (*this)(f(x[i]));

    return result;
}
gridfnc gridfnc::operator()(double (*f)(double))
{
    gridfnc result;

    result.set_gridpt(a, b, n);

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (*this)(f(x[i]));

    return result;
}

gridfnc& gridfnc::operator=(const gridfnc& f)
{
    x = f.x;
    y = f.y;
    a = f.a;
    b = f.b;
    h = f.h;
    n = f.n;

    return *this;
}
gridfnc& gridfnc::operator=(const double& val)
{
    y = val;
    return *this;
}
gridfnc& gridfnc::operator=(double (*f)(double))
{
    y = x.apply(f);
    return *this;
}

gridfnc& gridfnc::operator*= (const gridfnc& rhs)
{
    y *= rhs.y;
    return *this;
}
gridfnc& gridfnc::operator/= (const gridfnc& rhs)
{
    y /= rhs.y;
    return *this;
}
gridfnc& gridfnc::operator+= (const gridfnc& rhs)
{
    y += rhs.y;
    return *this;
}
gridfnc& gridfnc::operator-= (const gridfnc& rhs)
{
    y -= rhs.y;
    return *this;
}

gridfnc& gridfnc::operator*= (const double& val)
{
    y *= val;
    return *this;
}
gridfnc& gridfnc::operator/= (const double& val)
{
    y /= val;
    return *this;
}
gridfnc& gridfnc::operator+= (const double& val)
{
    y += val;
    return *this;
}
gridfnc& gridfnc::operator-= (const double& val)
{
    y -= val;
    return *this;
}

gridfnc& gridfnc::operator*= (double (*f)(double))
{
    y *= x.apply(f);
    return *this;
}
gridfnc& gridfnc::operator/= (double (*f)(double))
{
    y /= x.apply(f);
    return *this;
}
gridfnc& gridfnc::operator+= (double (*f)(double))
{
    y += x.apply(f);
    return *this;
}
gridfnc& gridfnc::operator-= (double (*f)(double))
{
    y -= x.apply(f);
    return *this;
}

// nonmember functions 
gridfnc operator* (const gridfnc& lhs, const gridfnc& rhs)
{
    gridfnc result = lhs;
    result.y *= rhs.y;

    return result;
}
gridfnc operator* (const double& val, const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y *= val;

    return result;
}
gridfnc operator* (const gridfnc& lhs, const double& val)
{
    gridfnc result = lhs;
    result.y *= val;

    return result;
}
gridfnc operator* (double (*f)(double), const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y *= result.x.apply(f);

    return result;
}
gridfnc operator* (const gridfnc& lhs, double (*f)(double))
{
    gridfnc result = lhs;
    result.y *= result.x.apply(f);

    return result;
}

gridfnc operator/ (const gridfnc& lhs, const gridfnc& rhs)
{
    gridfnc result = lhs;
    result.y /= rhs.y;

    return result;
}
gridfnc operator/ (const double& val, const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y = val/(rhs.y);

    return result;
}
gridfnc operator/ (const gridfnc& lhs, const double& val)
{
    gridfnc result = lhs;
    result.y /= val;

    return result;
}
gridfnc operator/ (double (*f)(double), const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y = (result.x.apply(f))/(rhs.y);

    return result;
}
gridfnc operator/ (const gridfnc& lhs, double (*f)(double))
{
    gridfnc result = lhs;
    result.y /= result.x.apply(f);

    return result;
}

gridfnc operator+ (const gridfnc& lhs, const gridfnc& rhs)
{
    gridfnc result = lhs;
    result.y += rhs.y;

    return result;
}
gridfnc operator+ (const double& val, const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y += val;

    return result;
}
gridfnc operator+ (const gridfnc& lhs, const double& val)
{
    gridfnc result = lhs;
    result.y += val;

    return result;
}
gridfnc operator+ (double (*f)(double), const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y += result.x.apply(f);

    return result;
}
gridfnc operator+ (const gridfnc& lhs, double (*f)(double))
{
    gridfnc result = lhs;
    result.y += result.x.apply(f);

    return result;
}

gridfnc operator- (const gridfnc& lhs, const gridfnc& rhs)
{
    gridfnc result = lhs;
    result.y -= rhs.y;

    return result;
}
gridfnc operator- (const double& val, const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y = val-(rhs.y);

    return result;
}
gridfnc operator- (const gridfnc& lhs, const double& val)
{
    gridfnc result = lhs;
    result.y -= val;

    return result;
}
gridfnc operator- (double (*f)(double), const gridfnc& rhs)
{
    gridfnc result = rhs;
    result.y = (result.x.apply(f)) - (rhs.y);

    return result;
}
gridfnc operator- (const gridfnc& lhs, double (*f)(double))
{
    gridfnc result = lhs;
    result.y -= result.x.apply(f);

    return result;
}

gridfnc max (const gridfnc& lhs, const gridfnc& rhs)
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] > rhs.y[i]) ? lhs.y[i] : rhs.y[i];

    return result;
}
gridfnc max (const double& val, const gridfnc& rhs)
{
    int n = rhs.n;
    gridfnc result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (val > rhs.y[i]) ? val : rhs.y[i];

    return result;
}
gridfnc max (const gridfnc& lhs, const double& val)
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] > val) ? lhs.y[i] : val;

    return result;
}
gridfnc max (double (*f)(double), const gridfnc& rhs)
{
    int n = rhs.n;
    gridfnc result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(rhs.x[i]) > rhs.y[i]) ? f(rhs.x[i]) : rhs.y[i];

    return result;
}
gridfnc max (const gridfnc& lhs, double (*f)(double))
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(lhs.x[i]) > lhs.y[i]) ? f(lhs.x[i]) : lhs.y[i];

    return result;
}

gridfnc min (const gridfnc& lhs, const gridfnc& rhs)
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] < rhs.y[i]) ? lhs.y[i] : rhs.y[i];

    return result;
}
gridfnc min (const double& val, const gridfnc& rhs)
{
    int n = rhs.n;
    gridfnc result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (val < rhs.y[i]) ? val : rhs.y[i];

    return result;
}
gridfnc min (const gridfnc& lhs, const double& val)
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] < val) ? lhs.y[i] : val;

    return result;
}
gridfnc min (double (*f)(double), const gridfnc& rhs)
{
    int n = rhs.n;
    gridfnc result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(rhs.x[i]) < rhs.y[i]) ? f(rhs.x[i]) : rhs.y[i];

    return result;
}
gridfnc min (const gridfnc& lhs, double (*f)(double))
{
    int n = lhs.n;
    gridfnc result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(lhs.x[i]) < lhs.y[i]) ? f(lhs.x[i]) : lhs.y[i];

    return result;
}

gridfnc convolve (const gridfnc& f, const gridfnc& g)
{
    gridfnc result = f;
    result.convolve_with( g );
    return result;
}
gridfnc convolve (const gridfnc& f, double (*g)(double))
{
    gridfnc result = f;
    result.convolve_with( g );
    return result;
}
gridfnc convolve (double (*g)(double), const gridfnc& f)
{
    gridfnc result = f;
    result.convolve_with( g );
    return result;
}

// overloaded cmath functions, 缺 atan2
gridfnc abs (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = abs(f.y);

    return result;
}
gridfnc acos (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = acos(f.y);

    return result;
}
gridfnc asin (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = asin(f.y);

    return result;
}
gridfnc atan (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = atan(f.y);

    return result;
}
gridfnc cos (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = cos(f.y);

    return result;
}
gridfnc cosh (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = cosh(f.y);

    return result;
}
gridfnc exp (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = exp(f.y);

    return result;
}
gridfnc log (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = log(f.y);

    return result;
}
gridfnc log10 (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = log10(f.y);

    return result;
}
gridfnc sin (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sin(f.y);

    return result;
}
gridfnc sinh (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sinh(f.y);

    return result;
}
gridfnc sqrt (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sqrt(f.y);

    return result;
}
gridfnc tan (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = tan(f.y);

    return result;
}
gridfnc tanh (const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = tanh(f.y);

    return result;
}
gridfnc pow (const gridfnc& f, const gridfnc& g)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(f.y, g.y);

    return result;
}
gridfnc pow (const gridfnc& f, const double& exponent)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(f.y, exponent);

    return result;
}
gridfnc pow (const double& base, const gridfnc& f)
{
    gridfnc result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(base, f.y);

    return result;
}

gridfnc ncdf (const gridfnc& f)
{
    // 為什麼不能直接 return f.apply( ncdf_double ); ?

    gridfnc result; 

    result.set_gridpt(f.a, f.b, f.n);
    result.y = f.y.apply( ncdf_double );

    return result;
}
gridfnc npdf (const gridfnc& f)
{
    gridfnc result; 

    result.set_gridpt(f.a, f.b, f.n);
    result.y = f.y.apply( npdf );

    return result;
}

#endif