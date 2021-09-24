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
//  template <class T> class gridfnc;
//  template <class T> class tablemaker
//  {
//  public:
//      void nbprint( const gridfnc<>& z );
//  };
//

using namespace std;

template <class T> class gridfnc
{
protected: 
    valarray<T> x, y;
    T a, b, h;
    int n;                                                            // number of grid points = n + 1
    static int order;
    static bool trap_conv;                                            // default value = 0; use midpt rule

    // for convolution
    valarray<T> repeatL (const valarray<T>& a);        // 在左邊抄 n-1 次的 a[0], n = a.size()
    valarray<T> repeatR (const valarray<T>& a);        // 在右邊抄 n-1 次的 a[n-1]，這兩個梯形法 convol 積分要用
    valarray<T> padR (const valarray<T>& a);        // 在右邊抄 n-1 次 0, for linear convolution
    valarray<T> truncate (const valarray<T>& a);    // 砍左右得到 convol 積分結果
    
public:
    gridfnc (){}
    gridfnc (const T& a, const T& b, const int& n, T (*f)(T));
    gridfnc (const T& a, const T& b, const int& n, const T& val);    // indicator 1_{x > val}
    gridfnc (const T& a, const T& b, const int& n);                    // f(x) = x 
    gridfnc (const gridfnc& f);                                        // copy constructor
    ~gridfnc (void){}
    
    static void set_order(int ord);
    static void use_trap_conv(bool flag=1);                            // use_trap_conv(): use trap rule for convolution
                                                                    // use_trap_conv(0): use midpt rule
    
    /*    ////////////////////////////////// 目前只有 2 階版本的成員函數 //////////////////////////////////
    T integral ();                                                // Trapezoidal rule
    T inner_prod (const gridfnc<T>& g);                            // integrate f*g from a to b using Trapezoidal rule, 假設 f.x = g.x
    T inner_prod (T (*g)(T));
    
    gridfnc<T> derivative ();                                    // 1st derivative (function), approximated by 2nd order finite difference
    T antiderivative ( const T& x );
    gridfnc<T> antiderivative ();

    void convolve_with (const gridfnc<T>& g);
    void convolve_with (T (*g)(T));
    
    ////////////////////////////////// 已經有 4 階版本的 //////////////////////////////////
    T operator() ( const T& val )
    T derivative ( const T& val )
    T derivative_2nd ( const T& val )
    */
    
    // set_gridpt 改變 x，a, b, n, h, y 也跟著變
    void set_gridpt (const T& a, const T& b, const int& n);            // set x grid pt, also used in constructors
    void set_gridpt (const gridfnc<T>& g);                            // set *this.x = g.x
    T integral ();                                                    // Trapezoidal rule
    T sup ();
    T inf ();
    T inner_prod (const gridfnc<T>& g);                                // integrate f*g from a to b using Trapezoidal rule, 假設 f.x = g.x
    T inner_prod (T (*g)(T));
    gridfnc<T> apply (T (*g)(T));                                    // g(*this) 反序合成，不需要 input gridfnc 的版本因為呼叫 () 就好了
    
    gridfnc<T> derivative ();                                        // 1st derivative (function), approximated by 2nd order finite difference
    T derivative ( const T& x );                                    // 1st derivative at x approximated by interpolating quadratic function
    T derivative_2nd ( const T& x );                                // 2nd derivative at x approximated by interpolating cubic function
    // gridfnc<T> antiderivative ();
    // T antiderivative ( const T& x );

    T root_near ( const T& x0 );                                    // Newton's method
    T root_between (T a, T b);                                        // bisection 
    T geta ();
    T getb ();
    T geth ();
    int getn ();
    int getindex(const T& val);                                        // find index i such that val is in (x[i], x[i+1])
    T getxi (const T& val);                                            // 回傳 val 左邊最近的格點
    T getyi (const T& val);                                            // 回傳 val 左邊最近的格點上的函數值
    valarray<T> getx ();
    valarray<T> gety ();
    void sety ( const valarray<T>& fncval );
    void setyi (int i, const T& val);

    // 假設兩函數用一樣的 x grid, 假設 a<0, b>0
    // 目前只能跑 T = double，如果要改 long double，要呼叫別的 fftw functions，如果要改 float，至少丟進去的 array 要先轉型成 double
    void convolve_with (const gridfnc<T>& g);
    void convolve_with (T (*g)(T));

    // operators
    gridfnc<T> operator+ () const;
    gridfnc<T> operator- () const;

    T operator() (const T& val);
    gridfnc<T> operator() (gridfnc<T>& f);            // 超載 () 時用 const 編譯不過。可以自己合成自己。g(h) 時用的是 g 的 x grid
    gridfnc<T> operator() (T (*f)(T));

    gridfnc<T>& operator= (const gridfnc<T>& f);
    gridfnc<T>& operator= (const T& val);
    gridfnc<T>& operator= (T (*f)(T));                // 不能丟內建函數，要自己另定義 double f(double) 才能用

    gridfnc<T>& operator*= (const gridfnc<T>& rhs);
    gridfnc<T>& operator/= (const gridfnc<T>& rhs);
    gridfnc<T>& operator+= (const gridfnc<T>& rhs);
    gridfnc<T>& operator-= (const gridfnc<T>& rhs);

    gridfnc<T>& operator*= (const T& val);
    gridfnc<T>& operator/= (const T& val);
    gridfnc<T>& operator+= (const T& val);
    gridfnc<T>& operator-= (const T& val);

    gridfnc<T>& operator*= (T (*f)(T));
    gridfnc<T>& operator/= (T (*f)(T));
    gridfnc<T>& operator+= (T (*f)(T));
    gridfnc<T>& operator-= (T (*f)(T));

    // nonmember functions                            operator+-*/ (T (*f)(T), const gridfnc<T>& rhs) 這種的四個都不能用
    //                                                神奇的是 max (T (*f)(T), const gridfnc<T>& rhs) 卻可以用
    //
    template <class T> friend gridfnc<T> operator* (const gridfnc<T>& lhs, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator* (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator* (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> operator* (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator* (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> operator/ (const gridfnc<T>& lhs, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator/ (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator/ (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> operator/ (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator/ (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> operator+ (const gridfnc<T>& lhs, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator+ (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator+ (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> operator+ (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator+ (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> operator- (const gridfnc<T>& lhs, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator- (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator- (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> operator- (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> operator- (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> max (const gridfnc<T>& lhs, const gridfnc<T>& rhs);    // 假設 f.x = g.x
    template <class T> friend gridfnc<T> max (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> max (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> max (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> max (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> min (const gridfnc<T>& lhs, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> min (const T& val, const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> min (const gridfnc<T>& lhs, const T& val);
    template <class T> friend gridfnc<T> min (T (*f)(T), const gridfnc<T>& rhs);
    template <class T> friend gridfnc<T> min (const gridfnc<T>& lhs, T (*f)(T));

    template <class T> friend gridfnc<T> convolve (const gridfnc<T>& f, const gridfnc<T>& g);    // 用 convolve_with，呼叫者會被改變；用 convolve 不會
    template <class T> friend gridfnc<T> convolve (const gridfnc<T>& f, T (*g)(T));
    template <class T> friend gridfnc<T> convolve (T (*g)(T), const gridfnc<T>& f);

    // overloaded cmath functions，缺 atan2
    template <class T> friend gridfnc<T> abs (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> acos (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> asin (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> atan (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> cos (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> cosh (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> exp (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> log (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> log10 (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> sin (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> sinh (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> sqrt (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> tan (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> tanh (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> pow (const gridfnc<T>& f, const gridfnc<T>& g);        // 假設 f.x = g.x
    template <class T> friend gridfnc<T> pow (const gridfnc<T>& f, const T& exponent);
    template <class T> friend gridfnc<T> pow (const T& base, const gridfnc<T>& f);

    // overloaded special functions
    template <class T> friend gridfnc<T> ncdf (const gridfnc<T>& f);
    template <class T> friend gridfnc<T> npdf (const gridfnc<T>& f);
};

template<class T> int gridfnc<T>::order = 2;
template<class T> bool gridfnc<T>::trap_conv = 0;

// private functions for convolution
template <class T> valarray<T> gridfnc<T>::repeatL(const valarray<T>& a)
{    
    int n = a.size();
    valarray <T> result(2*n-1);
    
    for(int i=0 ; i<n-1 ; i++)
        result[i] = a[0];
    for(int i=n-1 ; i<2*n-1 ; i++)
        result[i] = a[i-n+1];

    return result;
}
template <class T> valarray<T> gridfnc<T>::repeatR(const valarray<T>& a)
{
    int n = a.size();
    valarray <T> result(2*n-1);
    
    for(int i=0 ; i<n-1 ; i++)
        result[i] = a[i];
    for(int i=n-1 ; i<2*n-1 ; i++)
        result[i] = a[n-1];

    return result;
}
template <class T> valarray<T> gridfnc<T>::padR(const valarray<T>& a)
{
    int n = a.size();
    valarray<T> result(2*n-1);
    
    for( int i=0 ; i<a.size() ; i++ )
        result[i] = a[i];

    return result;
}
template <class T> valarray<T> gridfnc<T>::truncate(const valarray<T>& a)
{
    // a.size() is assumed to be odd
    // n is the size of the original array before padding 

    int n = (a.size() + 1)/2;
    valarray<T> result(n);
    
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
template <class T> gridfnc<T>::gridfnc(const T& a, const T& b, const int& n, T (*f)(T))
{
    set_gridpt(a, b, n);
    y = x.apply(f);
}
template <class T> gridfnc<T>::gridfnc(const T& a, const T& b, const int& n, const T& val)
{
    set_gridpt(a, b, n);
    // y = val; constant function 版本，後來改為 indicator 1_{x > val}

    for( int i=0 ; i<=n ; i++ )
        y[i] = (x[i]>val) ? 1 : 0;
}
template <class T> gridfnc<T>::gridfnc(const T& a, const T& b, const int& n)
{
    set_gridpt(a, b, n);
    y = x;
}
template <class T> gridfnc<T>::gridfnc(const gridfnc<T>& f)
{
    // f 的所有東西都是 private，可是拿的到
    
    set_gridpt(f.a, f.b, f.n);
    y = f.y;
}

template <class T> void gridfnc<T>::set_order(int ord)
{
    order = ord;
}
template <class T> void gridfnc<T>::use_trap_conv(bool flag=1)
{
    trap_conv = flag;
}

template <class T> void gridfnc<T>::set_gridpt (const T& lb, const T& ub, const int& num)
{
    a = lb;
    b = ub;
    n = num;
    h = (b-a)/n;
    
    if( x.size()!=0 )    // 本來就有 x, y 了。會呼叫 operator() 所以 x 和 y 是不能一邊算一邊蓋掉的
    {    
        valarray<T> tmpx(n+1);
        valarray<T> tmpy(n+1);
    
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
template <class T> void gridfnc<T>::set_gridpt (const gridfnc<T>& g)
{
    valarray<T> tmpy(g.n + 1);

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
template <class T> T gridfnc<T>::integral()
{
    return h*(y.sum() - y[0]/2 - y[n]/2);
}
template <class T> T gridfnc<T>::sup()
{
    // avoid discontinuity caused by numerical err at a, b
    valarray<T> tmpy = y;
    
    tmpy[0] = tmpy[1];
    tmpy[n] = tmpy[n-1];
    
    return tmpy.max();
}
template <class T> T gridfnc<T>::inf()
{
    // avoid discontinuity caused by numerical err at a, b
    valarray<T> tmpy = y;
    
    tmpy[0] = tmpy[1];
    tmpy[n] = tmpy[n-1];

    return tmpy.min();
}
template <class T> T gridfnc<T>::inner_prod(const gridfnc<T>& g)
{
    return h*(((g.y)*y).sum() - (g.y[0])*y[0]/2 - (g.y[n])*y[n]/2);
}
template <class T> T gridfnc<T>::inner_prod(T (*g)(T))
{
    return h*((x.apply(g)*y).sum() - (g(a))*y[0]/2 - (g(b))*y[n]/2);
}
template <class T> gridfnc<T> gridfnc<T>::apply(T (*g)(T))
{
    gridfnc<T> result; 

    result.set_gridpt(a, b, n);
    result.y = y.apply(g);

    return result;
}

template <class T> T gridfnc<T>::derivative( const T& val )
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

                T fpp[] = { (y[i-3] - 8*y[i-2] + 8*y[i]   - y[i+1])/(12*h), 
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
template <class T> T gridfnc<T>::derivative_2nd( const T& val )
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

                T fpp[] = {    (-y[i-3] + 16*y[i-2] -30*y[i-1] + 16*y[i]   - y[i+1])/(12*h*h), 
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
template <class T> gridfnc<T> gridfnc<T>::derivative()
{
    gridfnc<T> result;

    result.set_gridpt(a, b, n);

    for(int i=1 ; i<n ; i++)
        result.y[i] = (y[i+1] - y[i-1])/(2*h);

    result.y[0] = (-3*y[0] + 4*y[1] - y[2])/(2*h);
    result.y[n] = (3*y[n] - 4*y[n-1] + y[n-2])/(2*h);

    return result;
}

template <class T> T gridfnc<T>::root_near( const T& x0 )
{
    T tmpx, x = x0;
    T absTol = 0.000000000001;    // 10^{-12}
    T relTol = 0.000000000001;    // 10^{-12}

    //gridfnc<T> fp = derivative();

    for(int i=0 ; i<100 ; i++)
    {
        tmpx = x;            // need last x for increment test 
        // x = x - (*this)(x)/fp(x);
        T fp = derivative(x);
        x = x - (*this)(x)/fp;
        if( sign(abs(x-tmpx)-(absTol + relTol*abs(x))) < 0 )
            break;
    }

    return x;
}
template <class T> T gridfnc<T>::root_between (T a, T b)
{
    T x;
    T absTol = 0.000000000001;    // 10^{-12}
    T relTol = 0.000000000001;    // 10^{-12}
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

template <class T> int gridfnc<T>::getindex(const T& val)
{
    return (int)((val - x[0])/h);
}
template <class T> T gridfnc<T>::geta()
{
    return a;
}
template <class T> T gridfnc<T>::getb()
{
    return b;
}
template <class T> T gridfnc<T>::geth()
{
    return h;
}
template <class T> T gridfnc<T>::getxi(const T& val)
{
    return x[ getindex(val) ];
}
template <class T> T gridfnc<T>::getyi(const T& val)
{
    return y[ getindex(val) ];
}
template <class T> int gridfnc<T>::getn()
{
    return n;
}
template <class T> valarray<T> gridfnc<T>::getx ()
{
    return x;
}
template <class T> valarray<T> gridfnc<T>::gety ()
{
    return y;
}
template <class T> void gridfnc<T>::sety ( const valarray<T>& fncval )
{
    y = fncval;
}
template <class T> void gridfnc<T>::setyi ( int i, const T& val )
{
    y[i] = val;
}

// convolution
template <class T> void gridfnc<T>::convolve_with(const gridfnc<T>& g)
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
template <class T> void gridfnc<T>::convolve_with(T (*g)(T))
{
    gridfnc<T> tmp(a, b, n, g);
    (*this).convolve_with(tmp);
}

// operators
template <class T> gridfnc<T> gridfnc<T>::operator+() const
{
    return *this;
}
template <class T> gridfnc<T> gridfnc<T>::operator-() const
{
    gridfnc<T> result;

    result.set_gridpt(a, b, n);
    result.y = -y;
    return result;
}

template <class T> T gridfnc<T>::operator()(const T& val)
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
template <class T> gridfnc<T> gridfnc<T>::operator()(gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(a, b, n);

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (*this)(f(x[i]));

    return result;
}
template <class T> gridfnc<T> gridfnc<T>::operator()(T (*f)(T))
{
    gridfnc<T> result;

    result.set_gridpt(a, b, n);

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (*this)(f(x[i]));

    return result;
}

template <class T> gridfnc<T>& gridfnc<T>::operator=(const gridfnc<T>& f)
{
    x = f.x;
    y = f.y;
    a = f.a;
    b = f.b;
    h = f.h;
    n = f.n;

    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator=(const T& val)
{
    y = val;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator=(T (*f)(T))
{
    y = x.apply(f);
    return *this;
}

template <class T> gridfnc<T>& gridfnc<T>::operator*= (const gridfnc<T>& rhs)
{
    y *= rhs.y;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator/= (const gridfnc<T>& rhs)
{
    y /= rhs.y;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator+= (const gridfnc<T>& rhs)
{
    y += rhs.y;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator-= (const gridfnc<T>& rhs)
{
    y -= rhs.y;
    return *this;
}

template <class T> gridfnc<T>& gridfnc<T>::operator*= (const T& val)
{
    y *= val;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator/= (const T& val)
{
    y /= val;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator+= (const T& val)
{
    y += val;
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator-= (const T& val)
{
    y -= val;
    return *this;
}

template <class T> gridfnc<T>& gridfnc<T>::operator*= (T (*f)(T))
{
    y *= x.apply(f);
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator/= (T (*f)(T))
{
    y /= x.apply(f);
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator+= (T (*f)(T))
{
    y += x.apply(f);
    return *this;
}
template <class T> gridfnc<T>& gridfnc<T>::operator-= (T (*f)(T))
{
    y -= x.apply(f);
    return *this;
}

// nonmember functions 
template <class T> gridfnc<T> operator* (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    gridfnc<T> result = lhs;
    result.y *= rhs.y;

    return result;
}
template <class T> gridfnc<T> operator* (const T& val, const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y *= val;

    return result;
}
template <class T> gridfnc<T> operator* (const gridfnc<T>& lhs, const T& val)
{
    gridfnc<T> result = lhs;
    result.y *= val;

    return result;
}
template <class T> gridfnc<T> operator* (T (*f)(T), const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y *= result.x.apply(f);

    return result;
}
template <class T> gridfnc<T> operator* (const gridfnc<T>& lhs, T (*f)(T))
{
    gridfnc<T> result = lhs;
    result.y *= result.x.apply(f);

    return result;
}

template <class T> gridfnc<T> operator/ (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    gridfnc<T> result = lhs;
    result.y /= rhs.y;

    return result;
}
template <class T> gridfnc<T> operator/ (const T& val, const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y = val/(rhs.y);

    return result;
}
template <class T> gridfnc<T> operator/ (const gridfnc<T>& lhs, const T& val)
{
    gridfnc<T> result = lhs;
    result.y /= val;

    return result;
}
template <class T> gridfnc<T> operator/ (T (*f)(T), const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y = (result.x.apply(f))/(rhs.y);

    return result;
}
template <class T> gridfnc<T> operator/ (const gridfnc<T>& lhs, T (*f)(T))
{
    gridfnc<T> result = lhs;
    result.y /= result.x.apply(f);

    return result;
}

template <class T> gridfnc<T> operator+ (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    gridfnc<T> result = lhs;
    result.y += rhs.y;

    return result;
}
template <class T> gridfnc<T> operator+ (const T& val, const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y += val;

    return result;
}
template <class T> gridfnc<T> operator+ (const gridfnc<T>& lhs, const T& val)
{
    gridfnc<T> result = lhs;
    result.y += val;

    return result;
}
template <class T> gridfnc<T> operator+ (T (*f)(T), const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y += result.x.apply(f);

    return result;
}
template <class T> gridfnc<T> operator+ (const gridfnc<T>& lhs, T (*f)(T))
{
    gridfnc<T> result = lhs;
    result.y += result.x.apply(f);

    return result;
}

template <class T> gridfnc<T> operator- (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    gridfnc<T> result = lhs;
    result.y -= rhs.y;

    return result;
}
template <class T> gridfnc<T> operator- (const T& val, const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y = val-(rhs.y);

    return result;
}
template <class T> gridfnc<T> operator- (const gridfnc<T>& lhs, const T& val)
{
    gridfnc<T> result = lhs;
    result.y -= val;

    return result;
}
template <class T> gridfnc<T> operator- (T (*f)(T), const gridfnc<T>& rhs)
{
    gridfnc<T> result = rhs;
    result.y = (result.x.apply(f)) - (rhs.y);

    return result;
}
template <class T> gridfnc<T> operator- (const gridfnc<T>& lhs, T (*f)(T))
{
    gridfnc<T> result = lhs;
    result.y -= result.x.apply(f);

    return result;
}

template <class T> gridfnc<T> max (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] > rhs.y[i]) ? lhs.y[i] : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> max (const T& val, const gridfnc<T>& rhs)
{
    int n = rhs.n;
    gridfnc<T> result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (val > rhs.y[i]) ? val : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> max (const gridfnc<T>& lhs, const T& val)
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] > val) ? lhs.y[i] : val;

    return result;
}
template <class T> gridfnc<T> max (T (*f)(T), const gridfnc<T>& rhs)
{
    int n = rhs.n;
    gridfnc<T> result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(rhs.x[i]) > rhs.y[i]) ? f(rhs.x[i]) : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> max (const gridfnc<T>& lhs, T (*f)(T))
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(lhs.x[i]) > lhs.y[i]) ? f(lhs.x[i]) : lhs.y[i];

    return result;
}

template <class T> gridfnc<T> min (const gridfnc<T>& lhs, const gridfnc<T>& rhs)
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] < rhs.y[i]) ? lhs.y[i] : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> min (const T& val, const gridfnc<T>& rhs)
{
    int n = rhs.n;
    gridfnc<T> result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (val < rhs.y[i]) ? val : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> min (const gridfnc<T>& lhs, const T& val)
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (lhs.y[i] < val) ? lhs.y[i] : val;

    return result;
}
template <class T> gridfnc<T> min (T (*f)(T), const gridfnc<T>& rhs)
{
    int n = rhs.n;
    gridfnc<T> result;
    result.set_gridpt(rhs.a, rhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(rhs.x[i]) < rhs.y[i]) ? f(rhs.x[i]) : rhs.y[i];

    return result;
}
template <class T> gridfnc<T> min (const gridfnc<T>& lhs, T (*f)(T))
{
    int n = lhs.n;
    gridfnc<T> result;
    result.set_gridpt(lhs.a, lhs.b, n); 

    for(int i=0 ; i<=n ; i++)
        result.y[i] = (f(lhs.x[i]) < lhs.y[i]) ? f(lhs.x[i]) : lhs.y[i];

    return result;
}

template <class T> gridfnc<T> convolve (const gridfnc<T>& f, const gridfnc<T>& g)
{
    gridfnc<T> result = f;
    result.convolve_with( g );
    return result;
}
template <class T> gridfnc<T> convolve (const gridfnc<T>& f, T (*g)(T))
{
    gridfnc<T> result = f;
    result.convolve_with( g );
    return result;
}
template <class T> gridfnc<T> convolve (T (*g)(T), const gridfnc<T>& f)
{
    gridfnc<T> result = f;
    result.convolve_with( g );
    return result;
}

// overloaded cmath functions, 缺 atan2
template <class T> gridfnc<T> abs (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = abs(f.y);

    return result;
}
template <class T> gridfnc<T> acos (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = acos(f.y);

    return result;
}
template <class T> gridfnc<T> asin (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = asin(f.y);

    return result;
}
template <class T> gridfnc<T> atan (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = atan(f.y);

    return result;
}
template <class T> gridfnc<T> cos (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = cos(f.y);

    return result;
}
template <class T> gridfnc<T> cosh (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = cosh(f.y);

    return result;
}
template <class T> gridfnc<T> exp (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = exp(f.y);

    return result;
}
template <class T> gridfnc<T> log (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = log(f.y);

    return result;
}
template <class T> gridfnc<T> log10 (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = log10(f.y);

    return result;
}
template <class T> gridfnc<T> sin (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sin(f.y);

    return result;
}
template <class T> gridfnc<T> sinh (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sinh(f.y);

    return result;
}
template <class T> gridfnc<T> sqrt (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = sqrt(f.y);

    return result;
}
template <class T> gridfnc<T> tan (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = tan(f.y);

    return result;
}
template <class T> gridfnc<T> tanh (const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = tanh(f.y);

    return result;
}
template <class T> gridfnc<T> pow (const gridfnc<T>& f, const gridfnc<T>& g)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(f.y, g.y);

    return result;
}
template <class T> gridfnc<T> pow (const gridfnc<T>& f, const T& exponent)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(f.y, exponent);

    return result;
}
template <class T> gridfnc<T> pow (const T& base, const gridfnc<T>& f)
{
    gridfnc<T> result;

    result.set_gridpt(f.a, f.b, f.n);
    result.y = pow(base, f.y);

    return result;
}

template <class T> gridfnc<T> ncdf (const gridfnc<T>& f)
{
    // 為什麼不能直接 return f.apply( ncdf_double ); ?

    gridfnc<T> result; 

    result.set_gridpt(f.a, f.b, f.n);
    result.y = f.y.apply( ncdf_double );

    return result;
}
template <class T> gridfnc<T> npdf (const gridfnc<T>& f)
{
    gridfnc<T> result; 

    result.set_gridpt(f.a, f.b, f.n);
    result.y = f.y.apply( npdf );

    return result;
}

#endif