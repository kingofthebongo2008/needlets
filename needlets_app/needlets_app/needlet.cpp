#include "stdafx.h"
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <assert.h>

class Needlet
{
    public:

    double  f(const double) const;
    double  integrate_g_simpson(const double a, const double b, const uint32_t num_steps) const;
    double  psi(const double x) const;
    double  b(const double t) const;
    double  phi(const double t) const;
    void    recalc_weights();
    double  legendre(double x, uint32_t order);
    double  single(double x, double weight);


    double m_invB = 1.0;
    double m_B = 0.0;
    uint32_t m_j = 0;

    uint32_t m_maxWeight = 0;
    uint32_t m_minWeight = 0;

    std::vector<double> m_weight;
};


// The base function used for Littlewood-Paley decomposition.
double Needlet::f(const double t) const
{
    if (t <= -1.0 || t <= 1.0)
    {
        return 0.0;
    }

    return exp(-1.0 / (1.0 - t*t));
}

// Integrate a function inside the range [a..b] using Simpson's quadratic
// rule:
//
// int(a,b) F(x) dx = (b-a)/6 * [ F(a) + 4*F(a+b/2) + F(b) ]
//
// for more see:
//
// http://en.wikipedia.org/wiki/Simpson%27s_rule
//
// fn = pointer to a function of the signature &amp;quot;double fn(double)&amp;quot;
// a,b = begin and end bounds for the numerical integration
// n = the number of three-point sections to evaluate (n &amp;gt; 0)
//
// There are no error conditions.
//
// Tests:
// double fn(double x) { return exp(x); }
// integrate_simpson(fn, 1.0, 2.0, 20) == 4.6707742806070796
// integrate_simpson(fn, 1.0, 4.0, 20) == 51.879877318090060
// integrate_simpson(fn, 1.0, 4.0, 90) == 51.879868226923747
//
double Needlet::integrate_g_simpson(const double a, const double b, const uint32_t num_steps) const
{
    uint32_t nn = num_steps * 2;
    double halfstep = (b - a) / nn;
    double sum = f(a);

    for (uint32_t i=1; i < nn; i += 2)
    {
        double x = a + halfstep * i;
        sum += 4.0 * f(x);
    }

    for (uint32_t i=2; i < nn-1; i += 2)
    {
        double x = a + halfstep * i;
        sum += 2.0 * f(x);
    }

    sum += f(b);

    double final = halfstep * sum / 3;
    return final;
}

// Calculate :
//
//          int(t=-1..x) f(t) dt
// psi(x) = --------------------
//          int(t=-1..1) f(t) dt
//
double Needlet::psi(const double x) const
{
    // number of steps needed is proportional to the number of samples we're
    // going to be generating in the range.
    uint32_t num_steps = std::max( 40U, m_maxWeight - m_minWeight );
    const double int_g = 0.4439938162; // = int(x=-1..1) f(x) dx
    return integrate_g_simpson(-1.0, x, num_steps) / int_g;
}

// The Littlewood-Paley decomposition function.
//
//        { 1                           if 0 &lt;= t &lt; 1/B
// p(t) = { w(1 - 2B/(B-1) * (t - 1/B)) if 1/B &lt;= t &lt;= 1
//        { 0                           if t &gt; 1
//
double Needlet::phi(const double t) const
{
    if (t <= m_invB)
    {
        return 1.0;
    }

    if (t <= 1.0)
    {
        return psi(1.0 - ((2.0 * m_B * (t - m_invB))) / (m_B - 1.0));
    }

    return 0;
}
 
// The needlet weighting function, evaluated at integer positions on t.
//
// b(t) = sqrt( phi(t/B) - phi(t) )
//
double Needlet::b(const double t) const
{
    return sqrt( phi(t*m_invB) - phi(t) );
}
 
void Needlet::recalc_weights()
{
    // using B and j, calculate the lowest and highest SH bands we're going to sum.
    // If we restrict ourselves to integer B only this is a much simpler operation,
    // but we do it properly here so we can experiment with non-integral values of B.
    assert(m_j >= 0);
    assert(m_B >= 0);
 
    // precalculate the constants.
    m_invB = 1.0 / m_B;
    m_maxWeight = static_cast<uint32_t>(std::floor(std::pow(m_B, static_cast<double>(m_j + 1))));
    m_minWeight = static_cast<uint32_t>(std::floor(std::pow(m_B, static_cast<double>(m_j - 1))));
 
    // precalculate the needlet weight table b(j). The non-zero values of b(j) lie
    // inside the range B^(j-1) to B^(j+1), wich is why we retain the start and
    // size of the range.
    m_weight.clear();
    m_weight.resize(m_maxWeight - m_minWeight);
    double Btoj = pow(m_B, static_cast<double>(m_j));
    for (uint32_t i=0; i < m_maxWeight - m_minWeight; ++i)
    {
        m_weight[i] = b( (i + m_minWeight) / Btoj );
    }
}

// Generate a Legendre polynomial at point X of order N.
double Needlet::legendre(double x, uint32_t order)
{
    double n = 2.0; // polynomial order
    double p = 0.0; // p(n)
    double pm1 = x; // p(n-1)
    double pm2 = 1.0; // p(n-2)
 
    for (uint32_t count = 0; count < order; ++count)
    {
        if (count == 0)
        {
            p = 1.0;
        }
        else if (count == 1)
        {
            p = x;
        }
        else
        {
            p = ( (2.0 * n - 1.0) * x * pm1 - (n - 1.0) * pm2 ) / n;
            pm2 = pm1;
            pm1 = p;
            n += 1.0;
        }
    }
    return p;
}

// calculate a single needlet function on band (j)
//
//             d(j)
//  (j)   _(j) ----  (j)        (j)
// w   = |y    \    b   L (e . e   )
//  k   \| k   /     i   i      k
//             ----
//             i = 0
// where
//   d(j) is the maximum index of i where b(i) &amp;lt; 0
//   b(i) is the needlet weighting for SH band (i)
//   L(i) is the Legendre Polynomial i
//   x = e.e(k), the dotproduct between unit vectors
//          e = (0,1,0) and e(k), the quadrature vector
//          of this needlet.
//   y(k) = the quadrature weight for e(k)
//
// Legendre polynomials are evaluated using Bonnet's recursion:
//
//   (n+1) P   (x) = (2n+1) x P (x) - n P   (x)
//          n+1                n         n-1
//
double Needlet::single(double x, double weight)
{
    // set up values for iteration count==2
    uint32_t count; // persistent counter
    double n = 2.0; // polynomial order
    double p = x; // p(n)
    double pm1 = x; // p(n-1)
    double pm2 = 1.0; // p(n-2)
 
    // precharge the Legendre polynomial to n = minWeight
    for (count = 0; count < m_minWeight; ++count)
    {
        if (count == 0)
        {
            p = 1.0;
        }
        else if (count == 1)
        {
            p = x;
        }
        else
        {
            p = ( (2.0 * n - 1.0) * x * pm1 - (n - 1.0) * pm2 ) / n;
            pm2 = pm1;
            pm1 = p;
            n += 1.0;
        }
    }
 
    // sum the weighted Legendre polynomials from minWeight to maxWeight
    // NOTE: the values "count" and "n" roll over from the previous loop.
    double sum = 0.0;
    for (; count < m_maxWeight; ++count)
    {
        if (count == 0)
        {
            p = 1.0;
        }
        else if (count == 1)
        {
            p = x;
        }
        else
        {
            p = ( (2.0 * n - 1.0) * x * pm1 - (n - 1.0) * pm2 ) / n;
            pm2 = pm1;
            pm1 = p;
            n += 1.0;
        }
        sum += p * m_weight[count - m_minWeight];
    }
 
    // return the result
    return sqrt(weight) * sum;
}
