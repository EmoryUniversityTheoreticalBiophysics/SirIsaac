# plotMatrix.py
#
# Bryan Daniels
# 4.2.2012 moved from Sparseness.py
# 8.15.2013 moved from SparsenessTools.py
#

import scipy
import pylab
import mpl_toolkits.mplot3d.axes3d as p3 # for 3D plots

def plotMatrix(mat,cmap=pylab.cm.gray,colorbar=True,X=None,Y=None,          \
    interp='nearest',plot3D=False,plotContour=False,ax=None,                \
    contours=None,filled=True,outlineColor=None,**kwargs):
    """
    some popular cmaps:
        pylab.cm.gray
        pylab.cm.copper
        pylab.cm.bone
        pylab.cm.jet
        
    Can also use kwargs for pylab.imshow, including
        vmin,vmax: Set the range of values for the color bar.
        
    If 1D arrays for X and Y are given, axes are drawn.
    (unknown behavior if X and Y values are not equally spaced)
    Reminders:
        Y axis corresponds to 0 axis of matrix.
        X axis corresponds to 1 axis of matrix.
        Image will be flipped vertically compared to version w/o X and Y.
        
    For plot3D, useful kwargs:
        linewidth
    """
    if len(scipy.shape(mat)) == 1:
        mat = [mat]
    #minVal,maxVal = max(arrayFlatten(mat)),min(arrayFlatten(mat))
    if (ax is None) and plot3D:
      axNew = p3.Axes3D(pylab.figure())
    elif ax is None:
      pylab.figure()
      axNew = pylab.axes()
    else:
      pylab.axes(ax)
      axNew = ax
    
    if (plot3D or plotContour) and X is None: X = range(scipy.shape(mat)[1])
    if (plot3D or plotContour) and Y is None: Y = range(scipy.shape(mat)[0])
    
    if X is not None and Y is not None:
        X,Y = scipy.array(X),scipy.array(Y)
        
        # shift axes so they align with middle of each component
        if len(X) > 1:  deltaX = X[1] - X[0]
        else: deltaX = 1 
        if len(Y) > 1:  deltaY = Y[1] - Y[0]
        else: deltaY = 1
        if scipy.any(abs(X[1:]-X[:-1] - deltaX) > 1e-5) or                  \
           scipy.any(abs(Y[1:]-Y[:-1] - deltaY) > 1e-5):
           print "plotMatrix WARNING: X and/or Y values are not equally "+  \
                 "spaced.  May produce strange behavior."
          
        if plot3D or plotContour:
          
          Xmesh,Ymesh = scipy.meshgrid(X,Y)
          Z = scipy.array(mat)
          if plotContour:
            if filled:
              contourFn = axNew.contourf3D
            else:
              contourFn = axNew.contour3D
            if contours is None:
              contourLevels = contourFn(Xmesh,Ymesh,Z,extend3D=True,stride=1,**kwargs)
            else:
              contourLevels = contourFn(Xmesh,Ymesh,Z,contours,**kwargs)
          else:
            if filled:
              axNew.plot_surface(Xmesh,Ymesh,Z,rstride=1,cstride=1,         \
                cmap=cmap,**kwargs)
            else:
              axNew.plot_wireframe(Xmesh,Ymesh,Z,rstride=1,cstride=1,       \
                **kwargs)
        else:
          Xshifted = scipy.concatenate([[X[0]-deltaX],X]) + deltaX/2.
          Yshifted = scipy.concatenate([[Y[0]-deltaY],Y]) + deltaY/2.
          newAxis = [Xshifted[0],Xshifted[-1],Yshifted[0],Yshifted[-1]]
          if ax is not None:
            pylab.axes(ax)
            oldAxis = pylab.axis()
            newAxis = [min(newAxis[0],oldAxis[0]),                          \
                       max(newAxis[1],oldAxis[1]),                          \
                       min(newAxis[2],oldAxis[2]),                          \
                       max(newAxis[3],oldAxis[3])]             
          pylab.pcolor(Xshifted, Yshifted, mat, cmap=cmap, **kwargs)
          pylab.axis(newAxis)
          axNew = pylab.axes()
          if outlineColor is not None:
            rect = pylab.Rectangle([Xshifted[0],Yshifted[0]],
                            Xshifted[-1]-Xshifted[0],
                            Yshifted[-1]-Yshifted[0],fill=False,
                            ec=outlineColor,lw=0.25)
            axNew.add_patch(rect)
    else:
        #if ax is None:
        #  pylab.figure()
        #else:
        #  pylab.axes(ax)
        pylab.axes(axNew)
        pylab.imshow(mat, interpolation=interp, cmap=cmap, **kwargs)
        #axNew = pylab.axes()
        axNew.get_xaxis().set_visible(False)
        axNew.get_yaxis().set_visible(False)
        pylab.ylabel('')
        #ax = pylab.subplot(111)
        #ax.axis["right"].set_visible(False)
        #ax.axis["top"].set_visible(False)
        #grid(True)
    if colorbar and plotContour:
        pylab.colorbar(contourLevels,ax=axNew)
    elif colorbar:
        pylab.colorbar(ax=axNew)
        
    return axNew
    


