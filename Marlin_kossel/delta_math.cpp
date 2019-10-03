#include "Marlin.h"

//求3x3行列式
float determinant3x3(float n[3][3])
{
    return (
            n[0][0]*n[1][1]*n[2][2]+
            n[0][1]*n[1][2]*n[2][0]+
            n[0][2]*n[1][0]*n[2][1]-
            n[0][2]*n[1][1]*n[2][0]-
            n[0][0]*n[1][2]*n[2][1]-
            n[0][1]*n[1][0]*n[2][2]
            );
}
float d3x3(float n[3][3],float c[3],int cols)
{
    float s[3][3];
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            s[i][j]=n[i][j];
    for(int k=0;k<3;++k)
        s[k][cols]=c[k];
    return determinant3x3(s);
}
//交换a,b
void swap(float a[3],float b[3])
{
    float c[3];
    c[0]=b[0];
    c[1]=b[1];
    c[2]=b[2];
    b[0]=a[0];
    b[1]=a[1];
    b[2]=a[2];
    a[0]=c[0];
    a[1]=c[1];
    a[2]=c[2];
}
//解3元一次方程组
bool eq3x3(float n[3][3],float c[3],float x[3])
{
    float D = determinant3x3(n);
    if( D == 0 )
        return false; //无解
    else
    {
        for(int i=0;i<3;++i)
            x[i] = d3x3(n,c,i)/D;
        return true;
    }
}
/*
 通过点和法线构造平面方程
 平面上的一点p,和平面的法线n,构造平面方程a
 a[0]*x+a[1]*y+a[2]*z+a[3] = 0
 点法式方程
 n[0]*(x-p[0])+n[1]*(y-p[1])+n[2]*(z-p[2]) = 0
 */
void plane(float p[3],float n[3],float a[4])
{
    a[0] = n[0];
    a[1] = n[1];
    a[2] = n[2];
    a[3] = -(n[0]*p[0]+n[1]*p[1]+n[2]*p[2]);
}
//解一个元二次方程,将解放入x中，并且返回解的数量
int eq2(float a,float b,float c,float x[2])
{
    float e = sq(b)-4*a*c;
    if(e>=0)
    {
        if(e==0)
        {
            x[0] = -b/(2*a);
            return 1;
        }else
        {
            x[0] = (-b+sqrt(e))/(2*a);
            x[1] = (-b-sqrt(e))/(2*a);
            return 2;
        }
    }else return 0;
}
/*
 求解一个球面和直线的交点
 返回解的个数
 球面方程sphere[4],sq(x-sphere[0])+sq(y-sphere[1])+sq(z-sphere[2])=sphere[3]
 平面方程plane[4],plane[0]*x+plane[1]*y+plane[2]*z+plane[3] = 0
 球面方程s,平面方程p1,p2,返回的解c1,c2
 */
/*
    将一个两个平面组成的直线方程,转换为参数方程
    参数方程最后形式为:
    x = n[0]*t+p[0] ...
    方法使用平面z=0,z=1来求出两点,然后给出两点直线参数方程
    如果p1或者p2和xy平面平行,则使用x=0,x=1
 */
bool plane2line(float p1[4],float p2[4],float n[3],float p[3])
{
    float a[3][3],c[3],x[3];
    a[0][0] = p1[0];
    a[0][1] = p1[1];
    a[0][2] = p1[2];
    a[1][0] = p2[0];
    a[1][1] = p2[1];
    a[1][2] = p2[2];
    c[0] = -p1[3];
    c[1] = -p2[3];
    //增加平面z=0,然后求解
    a[2][0] = 0;
    a[2][1] = 0;
    a[2][2] = 1;
    c[2] = 0;
    if(eq3x3(a,c,p))
    {//如果有解
        //增加平面z=1,然后求解
        c[2] = 1;
        if(eq3x3(a,c,x))
        {
            n[0] = x[0]-p[0];
            n[1] = x[1]-p[1];
            n[2] = x[2]-p[2];
            return true;
        }
    }else
    { //方程无解,z=0和已知平面平行,该为x=0
        a[2][0] = 1;
        a[2][1] = 0;
        a[2][2] = 0;
        c[2] = 0;
        if(eq3x3(a, c, p))
        {
            c[2] = 1; //x=1
            if(eq3x3(a, c, x))
            {
                n[0] = x[0]-p[0];
                n[1] = x[1]-p[1];
                n[2] = x[2]-p[2];
                return true;
            }
        }
    }
    return false;
}
int sphere2plane(float s[4],float p1[4],float p2[4],float c1[3],float c2[3])
{
    float a,b,c,d,aa,bb,cc,x[2];
    aa = p1[1]*p2[2]-p2[1]*p1[2];
    if( aa==0 )return 0;
    a = (p2[0]*p1[2]-p1[0]*p2[2])/aa;
    b = (p2[3]*p1[2]-p1[3]*p2[2])/aa;
    aa = p1[2]*p2[1]-p2[2]*p1[1];
    if( aa==0 )return 0;
    c = (p2[0]*p1[1]-p1[0]*p2[1])/aa;
    d = (p2[3]*p1[1]-p1[3]*p2[1])/aa;
    
    aa = 1+sq(a)+sq(c);
    bb = 2*(-s[0]+a*(b-s[1])+c*(d-s[2]));
    cc = sq(s[0])+sq(b-s[1])+sq(d-s[2])-s[3];
    int ret = eq2(aa,bb,cc,x);
    if( ret>0 )
    {
        c1[0] = x[0];
        c1[1] = a*x[0]+b;
        c1[2] = c*x[0]+d;
        if( ret==2 )
        {
            c2[0] = x[1];
            c2[1] = a*x[1]+b;
            c2[2] = c*x[1]+d;
            return 2;
        }
        return 1;
    }
    return 0;
}
//方法2
int sphere2plane2(float s[4],float p1[4],float p2[4],float c1[3],float c2[3])
{
    float n[3],p[3],a,b,c,x[2];
    int ret;
    //首先转换为参数方程
    if(plane2line(p1, p2, n, p))
    {
        /*推导:直线x=n[0]*t+p[0],y=n[1]*t+p[1],z=n[2]*t+p[2]
         sq(x-s[0])+sq(y-s[1])+sq(z-s[2])=s[3]
         代入sq(n[0]*t+p[0]-s[0])+sq(n[1]*t+p[1]-s[1])+sq(n[2]*t+p[2]-s[2])=s[3]
         展开:(sq(n[0])+sq(n[1])+sq(n[2]))*sq(t) +
         2*(n[0]*(p[0]-s[0])+n[1]*(p[1]-s[1])+n[2]*(p[2]-s[2]))*t +
         sq(p[0]-s[0])+sq(p[1]-s[1])+sq(p[2]-s[2])-s[3] = 0
         */
        a = sq(n[0])+sq(n[1])+sq(n[2]);
        b = 2*(n[0]*(p[0]-s[0])+n[1]*(p[1]-s[1])+n[2]*(p[2]-s[2]));
        c = sq(p[0]-s[0])+sq(p[1]-s[1])+sq(p[2]-s[2])-s[3];
        ret = eq2(a,b,c,x);
        if( ret == 2 || ret == 1 )
        {
            for(int i=0;i<3;++i)
            {
                c1[i] = n[i]*x[0]+p[i];
                c2[i] = n[i]*x[1]+p[i];
            }
            return ret;
        }
    }
    return 0;
}
/*
 已知delta求点,calculate_delta的逆运算
 将求出的点放入destination
 */
bool calculate_cartesian(float delta[3],float cartesian[3])
{
    float o[3],d[3];
    float sqr,sqc,cosa,s;
    float sphere[4],p1[4],p2[4],c2[3];
    int c;
    o[X_AXIS] = (DELTA_TOWER1_X+DELTA_TOWER2_X)/2;
    o[Y_AXIS] = (DELTA_TOWER1_Y+DELTA_TOWER2_Y)/2;
    o[Z_AXIS] = (delta[X_AXIS]+delta[Y_AXIS])/2;
    
    sqr = sq(delta_diagonal_rod)-(sq(o[X_AXIS]-DELTA_TOWER1_X)+sq(o[Y_AXIS]-DELTA_TOWER1_Y)+sq(o[Z_AXIS]-delta[X_AXIS]));
    
    sqc = sq(o[X_AXIS]-DELTA_TOWER3_X)+sq(o[Y_AXIS]-DELTA_TOWER3_Y)+sq(o[Z_AXIS]-delta[Z_AXIS]);
    
    cosa = (sqc+sq(delta_diagonal_rod)-sqr)/(2*delta_diagonal_rod*sqrt(sqc));
    
    s = (cosa*delta_diagonal_rod)/sqrt(sqc);
    
    d[X_AXIS]=s*(o[X_AXIS]-DELTA_TOWER3_X)+DELTA_TOWER3_X;
    d[Y_AXIS]=s*(o[Y_AXIS]-DELTA_TOWER3_Y)+DELTA_TOWER3_Y;
    d[Z_AXIS]=s*(o[Z_AXIS]-delta[Z_AXIS])+delta[Z_AXIS];
    
    sphere[0] = DELTA_TOWER1_X;
    sphere[1] = DELTA_TOWER1_Y;
    sphere[2] = delta[X_AXIS];
    sphere[3] = sq(delta_diagonal_rod);
    
    c2[0] = DELTA_TOWER2_X-DELTA_TOWER1_X;
    c2[1] = DELTA_TOWER2_Y-DELTA_TOWER1_Y;
    c2[2] = delta[Y_AXIS]-delta[X_AXIS];
    plane(o,c2,p1);
    
    c2[0] = o[0]-DELTA_TOWER3_X;
    c2[1] = o[1]-DELTA_TOWER3_Y;
    c2[2] = o[2]-delta[Z_AXIS];
    plane(d,c2,p2);
    c=sphere2plane2(sphere,p1,p2,cartesian,c2);
    if( c>0 )
    {
        if(c==2&&cartesian[2]>c2[2])
        { //取z较小的解
            swap(cartesian,c2);
        }
        return true;
    }
    return false;
}
