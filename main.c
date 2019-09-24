#include <fcntl.h>
#include <math.h>
#include <omp.h>
#include <sys/param.h>

typedef struct Stack Stack;
typedef struct Queue Queue;
typedef struct Node Node;

struct Node {
    Node *next;
    unsigned int i;
    double d;
};

struct Queue {
    Node *first;
    Node *last;
};

struct Stack {
    Stack *head;
    unsigned int i;
};

int
put(Queue *q,
    const unsigned int i,
    const double d) {
    Node *n;
    n = malloc(sizeof(Node));
    if(!n)
        return 1;
    n->i = i;
    n->d = d;
    if(!q->first) {
        q->first = q->last = n;
    } else {
        q->last->next = n;
        q->last = n;
    }
    n->next = NULL;
    return 0;
}

int
get(Queue *q,
    unsigned int *i,
    double *d) {
    Node *tmp;
    if(!q->first)
        return 1;
    *i = q->first->i;
    *d = q->first->d;
    tmp = q->first;
    q->first = q->first->next;
    free(tmp);
    return 0;
}

unsigned int *
Simplicies(const unsigned int *tri,
           const unsigned int m,
           const unsigned int n) {
    unsigned int *net;
    unsigned int i, j, k, l, o;
    unsigned int *lst[n], ex;
    // node to simplicies map
    // m: number of facets
    // n: number of points

    // alloc list of arrays
    for(i = 0; i < n; i++) {
        lst[i] = malloc(8 * sizeof(unsigned int));
        if(!lst[i])
            return NULL;
        lst[i][0] = 8;
        lst[i][1] = 2;
    }

    l = 0;
    for(i = 0; i < m; i++) {
        for(j = 0; j < 3; j++) {
            k = tri[i*3+j];
            ex = 0;
            for(o = 2; o < lst[k][1]; o++) {
                if(lst[k][o] == i) {
                    ex = 1;
                    break;
                }
            }
            if(!ex) {
                lst[k][lst[k][1]++] = i;
                l++;
            }
            if(lst[k][0] < lst[k][1] + 2) {
                lst[k][0] = lst[k][1] + 4;
                lst[k] = realloc(lst[k], lst[k][0] * sizeof(unsigned int));
                if(!lst[k])
                    return NULL;
            }
        }
    }

    // store in compressed row format
    net = malloc((l + n + 1) * sizeof(unsigned int));
    if(!net)
        return NULL;
    j = n + 1;
    for(i = 0; i < n; i++) {
        net[i] = j;
        for(k = 2; k < lst[i][1]; k++) {
            net[j++] = lst[i][k];
        }
        free(lst[i]);
    }
    net[n] = j;
    return net;
}

unsigned int *
upstreamnetwork(const unsigned int *spx,
    const unsigned int m, unsigned int *r) {
    unsigned int *net;
    unsigned int i, j, k, l, o;
    unsigned int *lst[m], ex, itr;
    // reverse facet flow network spx

    // alloc list of arrays
    for(i = 0; i < m; i++) {
        lst[i] = malloc(8 * sizeof(unsigned int));
        if(!lst[i])
            return NULL;
        lst[i][0] = 8;
        lst[i][1] = 2;
    }

    l = 0;
    for(i = 0; i < m; i++) {
        itr = i * 2;
        for(j = 0; j < 2; j++) {
            k = spx[itr + j];
            if(k == m)
                continue;
            ex = 0;
            for(o = 2; o < lst[k][1]; o++) {
                if(lst[k][o] == i) {
                    ex = 1;
                    break;
                }
            }
            if(!ex) {
                lst[k][lst[k][1]++] = i;
                l++;
            }
            if(lst[k][0] < lst[k][1] + 2) {
                lst[k][0] = lst[k][1] + 4;
                lst[k] = realloc(lst[k], lst[k][0] * sizeof(unsigned int));
                if(!lst[k])
                    return NULL;
            }
        }
    }

    // store in compressed row format
    net = malloc((l + m + 1) * sizeof(unsigned int));
    if(!net)
        return NULL;
    j = m + 1;
    for(i = 0; i < m; i++) {
        net[i] = j;
        for(k = 2; k < lst[i][1]; k++) {
            net[j++] = lst[i][k];
        }
        free(lst[i]);
    }
    net[m] = j;
    *r = j;
    return net;
}

unsigned int
SimplexOfNodes(const unsigned int *net,
               const unsigned int a,
               const unsigned int b,
               const unsigned int x,
               const unsigned int m) {
    unsigned int i, k;

    // find the not-x simplex of two nodes a and b
    for(i = net[a]; i < net[a+1]; i++)
        for(k = net[b]; k < net[b+1]; k++)
            if(net[i] == net[k] && net[i] != x)
                return net[i];
    return m;
}

unsigned int
NodeOfSimplicies(const unsigned int *tri,
                 const unsigned int a,
                 const unsigned int b,
                 const double *z) {
    int r;
    unsigned int i, k, p, q;
    double zmin;

    // get the lowest node j of two facets a and b
    p = a*3;
    q = b*3;
    zmin = 1E99;
    r = -1;
    for(i = 0; i < 3; i++) {
        for(k = 0; k < 3; k++) {
            if(tri[p+i] == tri[q+k]) {
                if(z[tri[p+i]] < zmin) {
                    zmin = z[tri[p+i]];
                    r = tri[p+i];
                }
            }
        }
    }
    if(r == -1)
        exit(EXIT_FAILURE);
    return r;
}

double
HeronsTriangle(double a, double b, double c) {
    double d;
    // return the area of a facet

    // ! a >= b >= c
    if(a < b) {
        d = b;
        b = a;
        a = d;
    }
    if(b < c) {
        d = c;
        c = b;
        b = d;
    }
    if(a < b) {
        d = b;
        b = a;
        a = d;
    }
    return sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))) / 4;
}

void
linkthroughput(double *ltp,
               const unsigned int *spx,
               const double *spw,
               const double *spa,
               const unsigned int m) {
    double ltpi;
    unsigned int i, j, k, l;
    unsigned int *seen, *ideg, itr;
    Queue *que;

    // initialize
    seen = calloc(m, sizeof(unsigned int));
    ideg = calloc(m, sizeof(unsigned int));
    que = malloc(sizeof(Queue));
    if(!que || !ideg || !seen)
        exit(EXIT_FAILURE);
    que->first = que->last = NULL;

    // get in-degree
    for(i = 0; i < m; i++) {
        itr = i * 2;
        for(j = 0; j < 2; j++) {
            k = itr + j;
            l = spx[k];
            if(m > l)
                ideg[l]++;
        }
    }
    // start at facets without in-degree draining into l
    for(i = 0; i < m; i++) {
        if(!ideg[i]) {
            itr = i * 2;
            for(j = 0; j < 2; j++) {
                k = itr + j;
                l = spx[k];
                ltp[k] = spa[k];
                if(m > l)
                    if(put(que, l, ltp[k]))
                        exit(EXIT_FAILURE);
            }
        }
    }
    // work the queue
    while(!get(que, &i, &ltpi)) {
        seen[i]++;
        itr = i * 2;
        ltp[itr] += ltpi;
        if(seen[i] == ideg[i]) {
            // we collected all input for node i
            ltpi = ltp[itr];
            ltp[itr] = 0;
            for(j = 0; j < 2; j++) {
                k = itr + j;
                l = spx[k];
                // link throughput
                ltp[k] = ltpi * spw[k] + spa[k];
                if(m > l)
                    if(put(que, l, ltp[k]))
                        exit(EXIT_FAILURE);
            }
        }
    }
}

void
network(unsigned int *spx, double *spw, double *spa, double *spd,
	double *phi, double *theta,
        const unsigned int *tri,
        const double *x,
        const double *y,
        const double *z,
	const unsigned int m,
	const unsigned int n) {
    int sgn;
    double du, dv, dw, a, b, c;
    double xx, yy, slp, frc;
    double dx, dy, dz, dn, s, t;
    double xa, xb, xc, ya, yb, yc;
    double aa, ab, ac, bb, bc;
    double phii, beta;
    unsigned int i, j;
    unsigned int u, v, w, q, p;
    unsigned int *net;
    // m: number of facets
    // n: number of nodes

    net = Simplicies(tri, m, n);
    for(i = 0; i < m; i++) {
        // at p, q we store the pos of children
        p = i * 2;
        q = i * 2 + 1;
        for(j = 0; j < 3; j++) {
            u = tri[i*3 + j];
            v = tri[i*3 + (j+1)%3];
            w = tri[i*3 + (j+2)%3];
            // grad (dx,dy) of three point plane
            dz = ((x[w]-x[u])*(y[v]-y[u]) - (y[w]-y[u])*(x[v]-x[u]));
            dy = ((z[w]-z[u])*(x[v]-x[u]) - (x[w]-x[u])*(z[v]-z[u])) / dz;
            dx = ((y[w]-y[u])*(z[v]-z[u]) - (z[w]-z[u])*(y[v]-y[u])) / dz;

            // tri sides vs grad
            xa = x[w] - x[u];
            ya = y[w] - y[u];
            xb = x[v] - x[u];
            yb = y[v] - y[u];
            
            // dot products
            aa = xa*xa + ya*ya;
            ab = xa*xb + ya*yb;
            bb = xb*xb + yb*yb;
            dn = 1. / (aa*bb - ab*ab);
            for(sgn = -1; sgn <= 1; sgn += 2) {
                xc = sgn * dx;
                yc = sgn * dy;
                ac = xa*xc + ya*yc;
                bc = xb*xc + yb*yc;
                s = (bb*ac - ab*bc) * dn;
                t = (aa*bc - ab*ac) * dn;
                if(s >= 0 && t >= 0) {
                    phii = atan2(dy, dx);
                    phi[i] = phii;
		    theta[i] = atan(sqrt(dx*dx + dy*dy));
                    if(phii < 0)
                        phii += M_PI;
                    a = sqrt(xa*xa + ya*ya);
                    b = sqrt(xb*xb + yb*yb);
                    if(sgn > 0) {
                        spx[p] = SimplexOfNodes(net, w, v, i, m);
                        spx[q] = m;
                        spw[p] = 1;
                        spw[q] = 0;
                        c = sqrt((x[v]-x[w])*(x[v]-x[w])+(y[v]-y[w])*(y[v]-y[w]));
                        spa[p] = HeronsTriangle(a, b, c);
                        spa[q] = 0;
                        beta = atan2(y[w]-y[v], x[w]-x[v]);
                        if(beta < 0)
                            beta += M_PI;
                        beta -= phii;
                        if(beta > M_PI / 2)
                            beta = M_PI - beta;
                        spd[i] = c * fabs(sin(beta));
                    } else {
                        slp = dy / dx;
                        frc = (y[w] - y[v]) / (x[w] - x[v]);
                        if(dx) {
                            if(x[w] != x[v])
                                xx = (yb + x[u]*slp - x[v]*frc) / (slp - frc);
                            else
                                xx = x[w];
                            yy = (xx - x[u])*slp + y[u];
                        } else {
                            xx = x[u];
                            yy = (xx - x[w])*frc + y[w];
                        }
                        if(isinf(yy)) {
                            fprintf(stderr, "flat triangle %i (u:%.2f v:%.2f w:%.2f)\n", i, z[u], z[v], z[w]);
                            spw[p] = 0.5;
                            spw[q] = 0.5;
                            c = sqrt((x[v]-x[w])*(x[v]-x[w])+(y[v]-y[w])*(y[v]-y[w]));
                            spa[p] = HeronsTriangle(a, b, c) / 2.0;
                            spa[q] = spa[p];
                        } else {
                            du = sqrt((xx-x[u])*(xx-x[u])+(yy-y[u])*(yy-y[u]));
                            dv = sqrt((xx-x[v])*(xx-x[v])+(yy-y[v])*(yy-y[v]));
                            dw = sqrt((xx-x[w])*(xx-x[w])+(yy-y[w])*(yy-y[w]));
                            spw[p] = dv / (dv+dw);
                            spw[q] = dw / (dv+dw);
                            spa[p] = HeronsTriangle(b, dv, du);
                            spa[q] = HeronsTriangle(a, dw, du);
                        }
                        spx[p] = SimplexOfNodes(net, u, v, i, m);
                        spx[q] = SimplexOfNodes(net, u, w, i, m);
                        beta = atan2(yb, xb);
                        if(beta < 0)
                            beta += M_PI;
                        beta -= phii;
                        if(beta > M_PI / 2)
                            beta = M_PI - beta;
                        spd[i] = b * fabs(sin(beta));
                        beta = atan2(ya, xa);
                        if(beta < 0)
                            beta += M_PI;
                        beta -= phii;
                        if(beta > M_PI / 2)
                            beta = M_PI - beta;
                        spd[i] += a * fabs(sin(beta));
                    }
                    j = 3;
                    break;
                }
            }
        }
    }
    free(net);
}

void
tunnel(unsigned int *spx, double *spw, double *spa,
          const unsigned int *tri,
          const double *x,
          const double *y,
          const double *z,
          const unsigned int m,
          const unsigned int n,
          const double tubemaxdist) {
    double zu, zv, dv;
    unsigned int p, q, u, v, w;
    unsigned int i, j, k, l, s, t;
    unsigned int msinks, nsinks, mm;
    unsigned int *net, *seen, *sinks, dst;
    unsigned int *sinku, *uniqu, *udest;
    Queue *que;

    // m: number of facets
    // n: number of points
    net = Simplicies(tri, m, n);
    mm = m + m;

    sinks = malloc(mm * 2 * sizeof(unsigned int));
    sinku = malloc(mm * 2 * sizeof(unsigned int));
    if(!sinks || !sinku)
        exit(EXIT_FAILURE);

//#pragma omp parallel for private(i,j,k,l,s,p,q,u,v,zu,zv,dst)
//we don't want to have this in parallel because we manipulate spx[l*2+k] = dst
    for(i = 0; i < m; i++) {
        p = i * 2;
        for(j = 0; j < 2; j++) {
            q = p + j;
            sinks[q] = mm;
            sinku[q] = n;
            l = spx[q];
            if(l == m)
                continue;
            // check whether two neighboring facets flow into each other
            if(spx[l*2] == i || spx[l*2+1] == i) {
                // get lowest node of these two facets
                u = NodeOfSimplicies(tri, i, l, z);
                zu = z[u];
                dst = m;
                for(k = net[u]; k < net[u+1]; k++) {
                    v = net[k];
                    if(v == i || v == l)
                        continue;
                    zv = z[tri[v*3]];
                    if(z[tri[v*3+1]] > zv)
                        zv = z[tri[v*3+1]];
                    if(z[tri[v*3+2]] > zv)
                        zv = z[tri[v*3+2]];
                    if(zv == zu) {
                        dst = v;
                        break;
                    }
                }
                if(dst < m) {
                    spx[q] = dst;
                    // rewire also the other facet to that lower facet (l->dest)
                    for(k = 0; k < 2; k++)
                        if(spx[l*2+k] == i)
                            spx[l*2+k] = dst;
                } else {
                    sinks[q] = q;
                    sinku[q] = u;
                }
            }
        }
    }

    msinks = 0;
    for(i = 0; i < mm; i++) {
        if(sinks[i] < mm) {
            sinks[msinks] = sinks[i];
            sinku[msinks++] = sinku[i];
        }
    }
    sinks = realloc(sinks, msinks * sizeof(unsigned int));
    sinku = realloc(sinku, msinks * sizeof(unsigned int));
    uniqu = malloc(n * sizeof(unsigned int));
    udest = malloc(n * sizeof(unsigned int));
    if(!uniqu || !udest)
        exit(EXIT_FAILURE);

#pragma omp parallel for
    for(i = 0; i < n; i++) {
        uniqu[i] = n;
        udest[i] = m;
    }

    for(i = 0; i < msinks; i++) {
        u = sinku[i];
        uniqu[u] = u;
    }
    
    nsinks = 0;
    for(i = 0; i < n; i++)
        if(uniqu[i] < n)
            uniqu[nsinks++] = uniqu[i];
    uniqu = realloc(uniqu, nsinks * sizeof(unsigned int));
    
#pragma omp parallel for private(i,k,s,t,u,v,w,zu,zv,dv,seen,que) schedule(dynamic, 4)
    for(i = 0; i < nsinks; i++) {
        u = uniqu[i];
        zu = z[u];
        que = malloc(sizeof(Queue));
        seen = calloc(m, sizeof(unsigned int));
        if(!que || !seen)
            exit(EXIT_FAILURE);
        que->first = que->last = NULL;
        for(k = net[u]; k < net[u+1]; k++) {
            v = net[k];
            seen[v]++;
            if(put(que, v, 0))
                exit(EXIT_FAILURE);
        }
        while(!get(que, &v, &dv)) {
            if(dv > tubemaxdist) {
                break;
            }
            zv = z[tri[v*3]];
            if(z[tri[v*3+1]] > zv)
                zv = z[tri[v*3+1]];
            if(z[tri[v*3+2]] > zv)
                zv = z[tri[v*3+2]];
            if(zv < zu) {
                udest[u] = v;
                break;
            }
            for(s = 0; s < 3; s++) {
                t = tri[v*3+s];
                for(k = net[t]; k < net[t+1]; k++) {
                    w = net[k];
                    if(seen[w])
                        continue;
                    seen[w]++;
                    if(put(que, w, dv+1))
                        exit(EXIT_FAILURE);
                }
            }
        }
        while(!get(que, &v, &dv));
        free(seen);
        free(que);
    }
    free(uniqu);
    free(net);

    //tend = omp_get_wtime();
    //printf("%.4f\n", tend - tini);
    //tini = tend;

#pragma omp parallel for private(i,q,u)
    for(i = 0; i < msinks; i++) {
        q = sinks[i];
        u = sinku[i];
        spx[q] = udest[u];
    }
    free(udest);
    free(sinks);
    free(sinku);

    // fix spw and spa
#pragma omp parallel for private(i,p)
    for(i = 0; i < m; i++) {
        p = i * 2;
        if(spx[p] == m && spx[p+1] < m) {
            spa[p+1] += spa[p];
            spw[p+1] = 1;
            spw[p] = 0;
        } else if(spx[p] < m && spx[p+1] == m) {
            spa[p] += spa[p+1];
            spw[p] = 1;
            spw[p+1] = 0;
        }
    }
}

void
rivers(unsigned int *ind, const double *sca,
       const unsigned int *net, const unsigned int *rev,
       const unsigned int m, const double fac) {
    unsigned int i, j, k;
    Queue *que;
    double nl;

    que = malloc(sizeof(Queue));
    if(!que)
	exit(EXIT_FAILURE);
    que->first = que->last = NULL;
    for(i = 0; i < m; i++) {
	if(ind[i])
	    if(put(que, i, 0))
		exit(EXIT_FAILURE);
    }
    while(!get(que, &i, &nl)) {
	for(j = rev[i]; j < rev[i+1]; j++) {
	    k = rev[j];
	    if(ind[k])
		continue;
	    if(sca[k] < sca[i]*fac)
		continue;
	    ind[k] = 1;
	    if(put(que, k, 0))
                exit(EXIT_FAILURE);
	}
    }
    free(que);
}

void
convergence(double *conv, const double *sca,
            const unsigned int *net, const unsigned int *rev,
	    const unsigned int *sub, const unsigned int slen,
            const unsigned int m, const unsigned int nsamples) {
    unsigned int i, j, k, l, s, t, d;
    double qoff, qdia, sum, wgh, ak, aj, nl;
    Queue *que;
    unsigned int *dtr, *utr, *seen, *mask;

    mask = calloc(m, sizeof(unsigned int));
    if(!mask)
	exit(EXIT_FAILURE);
#pragma omp parallel for private(s)
    for(s = 0; s < slen; s++)
	mask[sub[s]] = 1;

#pragma omp parallel for private(i)
    for(i = 0; i < m; i++) {
	if(mask[i])
	    mask[i] = 0;
	else
	    mask[i] = 1;
    }

#pragma omp parallel for private(i,j,k,l,s,t,d,ak,aj,qoff,qdia,sum,wgh,nl,que,dtr,utr,seen)
    for(s = 0; s < slen; s++) {
	// downstream window
	i = sub[s];
	seen = malloc(m * sizeof(unsigned int));
	que = malloc(sizeof(Queue));
	dtr = malloc(nsamples * sizeof(unsigned int));
	if(!que || !seen || !dtr)
	    exit(EXIT_FAILURE);
	memcpy(seen, mask, m * sizeof(unsigned int));
        que->first = que->last = NULL;
	if(put(que, i, 0))
	    exit(EXIT_FAILURE);
	seen[i] = 1;
	d = 0;
	while(!get(que, &j, &nl)) {
            for(l = 0; l < 2; l++) {
                k = net[l+j*2];
		if(k == m)
		    continue;
		if(seen[k])
		    continue;
		seen[k] = 1;
                if(put(que, k, 0))
                    exit(EXIT_FAILURE);
		dtr[d++] = k;
		if(d == nsamples) {
		    while(!get(que, &j, &nl));
		    break;
		}
	    }
	}
	free(que);
	free(seen);
	if(d < nsamples) {
	    conv[i] = NAN;
	    free(dtr);
	    continue;
	}
	// upstream window
	seen = malloc(m * sizeof(unsigned int));
	que = malloc(sizeof(Queue));
	utr = malloc(nsamples * sizeof(unsigned int));
	if(!que || !seen || !utr)
	    exit(EXIT_FAILURE);
	memcpy(seen, mask, m * sizeof(unsigned int));
        que->first = que->last = NULL;
	if(put(que, i, 0))
	    exit(EXIT_FAILURE);
	seen[i] = 1;
	d = 0;
	while(!get(que, &j, &nl)) {
            for(l = rev[j]; l < rev[j+1]; l++) {
                k = rev[l];
		if(seen[k])
		    continue;
		seen[k] = 1;
                if(put(que, k, 0))
                    exit(EXIT_FAILURE);
		utr[d++] = k;
		if(d == nsamples) {
		    while(!get(que, &j, &nl));
		    break;
		}
	    }
	}
	free(que);
	free(seen);
	if(d < nsamples) {
	    conv[i] = NAN;
	    free(dtr);
	    free(utr);
	    continue;
	}
	sum = 0;
	wgh = 0;
	for(l = 1; l < d; l++) {
	    ak = sca[dtr[l]];
	    for(t = 0; t < l; t++) {
		aj = sca[dtr[t]];
		sum += abs(ak - aj);
		wgh += 1;
	    }
	}
	for(l = 1; l < d; l++) {
	    ak = sca[utr[l]];
	    for(t = 0; t < l; t++) {
		aj = sca[utr[t]];
		sum += abs(ak - aj);
		wgh += 1;
	    }
	}
	qdia = sum / wgh;
	sum = 0;
	wgh = 0;
	for(l = 0; l < d; l++) {
	    ak = sca[utr[l]];
	    for(t = 0; t < d; t++) {
		aj = sca[dtr[t]];
		sum += abs(ak - aj);
		wgh += 1;
	    }
	}
	qoff = sum / wgh;
	conv[i] = qoff - qdia;
	free(dtr);
	free(utr);
    }
}
