#!/usr/bin/env python


import matplotlib.pyplot as plt
import math


"""
Cubic Spline Interpolation
Generate a smooth (interpolated) curve that follows the path of the given X/Y points
https://swharden.com/blog/2022-01-22-spline-interpolation/
and
https://swharden.com/blog/2022-06-23-resample-interpolation/
"""


class SPLINE:

    def csvtolists(self,filename):
        x_points = []
        y_points = []
        f=open(filename,'r')
        first_row = True
        counter = 0
        for rows in f:
            # first row contains no numbers but the keys
            if first_row:
                first_row = False
                continue #stop current iteration, continue with next iteration
            # split CSV-rows
            numbers = rows.split(',')  # CSV
            #print(numbers)
            try:
                x = float(numbers[0])
            except:
                x = 0
            try:
                y = float(numbers[1])
            except:
                y = 0
            x_points.append(x)
            y_points.append(y)
            counter += 1
        print("List with %d datapoints generated." %counter)
        return (x_points, y_points)

        
    def FitMatrix(self, x, y):
        n = len(x)
        a = [0.0] * (n-1)
        b = [0.0] * (n-1)
        r = [0.0] * n
        A = [0.0] * n
        B = [0.0] * n
        C = [0.0] * n
                
        dx1 = x[1] -x[0]
        C[0] = 1.0 / dx1
        B[0] = 2.0 * C[0]
        r[0] = 3 * (y[1] - y[0])/(dx1 * dx1)

        for i in range(1, n-1, 1):
            dx1 = x[i] - x[i - 1]
            dx2 = x[i + 1] - x[i]
            A[i] = 1.0 / dx1
            C[i] = 1.0 / dx2
            B[i] = 2.0 * (A[i] + C[i])
            dy1 = y[i] - y[i - 1]
            dy2 = y[i + 1] - y[i]
            r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2))
        
        dx1 = x[n - 1] - x[n - 2]
        dy1 = y[n - 1] - y[n - 2]
        A[n - 1] = 1.0 / dx1
        B[n - 1] = 2.0 * A[n - 1]
        r[n - 1] = 3 * (dy1 / (dx1 * dx1))
        #print(x)
        #print(y)

        cPrime = [0.0] * n
        cPrime[0] = C[0] / B[0]
        for i in range(1, n, 1):
            cPrime[i] = C[i] / (B[i] - cPrime[i - 1] * A[i])
        
        dPrime = [0.0] * n
        dPrime[0] = r[0] / B[0]
        for i in range(1, n, 1):
            dPrime[i] = (r[i] - dPrime[i - 1] * A[i]) / (B[i] - cPrime[i - 1] * A[i])
        
        k = [0.0] * n
        k[n - 1] = dPrime[n - 1]
        for i in range(n-2, 0, -1):
            k[i] = dPrime[i] - cPrime[i] * k[i + 1]
        
        for i in range(1, n, 1):    
            dx1 = x[i] - x[i - 1]
            dy1 = y[i] - y[i - 1]
            a[i - 1] = k[i - 1] * dx1 - dy1
            b[i - 1] = -k[i] * dx1 + dy1
        
        return(a,b)
    

    def Interpolate(self, xOrig, yOrig, xInterp):
        (a,b) = SPLINE.FitMatrix(self, xOrig, yOrig)
        yInterp = [None] * len(xInterp)

        for i in range(0, len(yInterp), 1):
            for j in range(0, len(xInterp) -2, 1):
                if (xInterp[i] <= xOrig[j + 1]):
                    break

            dx = xOrig[j + 1] - xOrig[j]
            t = (xInterp[i] - xOrig[j]) / dx
            y = (1 - t) * yOrig[j] + t * yOrig[j + 1] + t * (1 - t) * (a[j] * (1 - t) + b[j] * t)
            yInterp[i] = y
        
        return yInterp
    

    def InterpolateXY(self, xs, ys, count):
        if (xs == 0 or ys == 0 or len(xs) != len(ys)):
            print("xs and ys must have same length")
        
        inputPointCount = len(xs)
        inputDistances = [0.0] * inputPointCount
        
        for i in range(1, inputPointCount, 1):
            dx = xs[i] - xs[i - 1]
            dy = ys[i] - ys[i - 1]
            distance = math.sqrt(dx * dx + dy * dy)
            inputDistances[i] = inputDistances[i - 1] + distance

        #print("inputDistances: ", inputDistances)
        
        meanDistance = inputDistances[-1] / (count - 1)
        #print("meanDistance: ", meanDistance)
        #evenDistances = Enumerable.Range(0, count).Select(x => x * meanDistance).ToArray();
        evenDistances = [0.0] * count
        for i in range(0, count, 1):
            evenDistances[i] = i * meanDistance
        
        #print("evenDistances: ", evenDistances)

        xsOut = SPLINE.Interpolate(self, inputDistances, xs, evenDistances)
        ysOut = SPLINE.Interpolate(self, inputDistances, ys, evenDistances)
        return (xsOut, ysOut)
    

    def Interpolate1D(self, xs, ys, count):
        if (xs == 0 or ys == 0 or len(xs) != len(ys)):
            print("xs and ys must have same length")
        
        inputPointCount = len(xs)
        inputDistances = [0.0] * inputPointCount
        
        for i in range(1, inputPointCount, 1):
            inputDistances[i] = inputDistances[i - 1] + xs[i] - xs[i - 1]

        #print("inputDistances: ", inputDistances)
        
        # letztes x-Element / (Anzahl Strecken = Anzahl Stützpunkte - 1)
        meanDistance = inputDistances[-1] / (count - 1)
        #print("meanDistance: ", meanDistance)
        #evenDistances = Enumerable.Range(0, count).Select(x => x * meanDistance).ToArray();
        evenDistances = [0.0] * count
        for i in range(0, count, 1):
            evenDistances[i] = i * meanDistance

        
        #print("evenDistances: ", evenDistances)

        xsOut = SPLINE.Interpolate(self, inputDistances, xs, evenDistances)
        ysOut = SPLINE.Interpolate(self, inputDistances, ys, evenDistances)
        return (xsOut, ysOut)


    # plot Array in diagram
    def plot_array(self, x_sample, y_sample, xpoints, ypoints, xpoints2, ypoints2, filename):
        # Plot
        plt.figure(figsize=(8, 6))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.plot(x_sample, y_sample, 'o', label='Sample', markersize=10.0)
        plt.plot(xpoints, ypoints, label='Cubic Spline', marker="o", markersize=5.0)
        plt.plot(xpoints2, ypoints2, label='Original', marker="x", markersize=5.0)
        plt.legend(['Sample', 'Cubic Spline', "Original"])
        plt.title(filename)
        plt.margins(0, .05)
        plt.tight_layout()
        plt.grid()
        plt.savefig(filename) # save Plot
        plt.show()


    # Sinusfunktion
    def sinusfct(self):
        x_list =[0.0] * 361
        y_list =[0.0] * 361
        for i in range(0,361,1):
            x_list[i] = i
            y_list[i] = math.sin(i*math.pi/180)
        #print(y_list)
        return(x_list, y_list)
        

def main():
    werte=SPLINE()
    x_sample = [0,45,135,180,225,270,315,360]
    y_sample = [0,0.707,0.707,0,-0.707,-1,-0.707,0]
    
    x_sample = [0,44.43,133,180,223.1,269.5,311.41,360]
    y_sample = [0,0.7,0.731,0,-0.6833,-0.999,-0.75,0]
    
    filename = "SplineTest.csv"
    (x_sample, y_sample) = werte.csvtolists(filename)
    print(x_sample,y_sample)

    (x_orig, y_orig) = werte.sinusfct()

    (x,y) = werte.Interpolate1D(x_sample, y_sample, 3601)
    #(x,y) = werte.InterpolateXY(x_sample, y_sample, 361)


    werte.plot_array(x_sample, y_sample, x, y, x_orig, y_orig, "SplineTest")

    
    
if __name__ == "__main__":
    main()