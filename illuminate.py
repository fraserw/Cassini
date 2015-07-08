import numpy as num

def midpoint(a,b,c):
    return num.array([num.mean([a[0],b[0],c[0]]),
                num.mean([a[1],b[1],c[1]]),
                num.mean([a[2],b[2],c[2]])])

def mag(v):
    return num.sum(v**2)**0.5

def getNormalMidpoint(p0,p1,p2):
    m=midpoint(p0,p1,p2)
    #n_centre_midpoint=m/mag(m)

    a=p0-m
    a/=mag(a)
    b=p1-m
    b/=mag(b)
    c=p2-m
    c/=mag(c)

    
    #wihtout negatives, seems to make normals that point outwards
    n_ab=num.cross(a,b)
    n_bc=num.cross(b,c)
    n_ca=num.cross(c,a)

    #seems to make normals that point inwards
    #n_ab=num.cross(b,a)
    #n_bc=num.cross(c,b)
    #n_ca=num.cross(a,c)
    
    n=midpoint(n_ab/mag(n_ab),n_bc/mag(n_bc),n_ca/mag(n_ca))
    print n
    #make simplification to take -n if n-m if further away from origin
    #than n+m
    x=m-n
    nx=m+n
    if mag(x)>mag(nx): n*=-1

    
    #faceAngle_dot=num.dot(n,n_centre_midpoint)
    #if faceAngle_dot<0:
    #    print faceAngle_dot,
    
    #if mag(n_ab)==0. or mag(n_bc)==0. or mag(n_ca)==0.:
    #    print n
    return (n,m)


