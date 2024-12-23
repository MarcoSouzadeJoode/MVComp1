#
# Simple viewing tool for arrays
#
# Copyright (c) 2018 C.P. Dullemond
# Free software under the standard MIT License
#
def viewarr(data,index=0,x=None,ymin=None,ymax=None,ylabel=None,idxnames=None,idxvals=None,idxformat=''):
    """
    Interactive plot of a 1-D cut from an n-dimensional array.

    EXAMPLE 1:
    from viewarr import *
    data=np.arange(64).reshape((4,4,4)) # Dummy dataset
    viewarr(data)

    EXAMPLE 2:
    from viewarr import *
    data=np.arange(64).reshape((4,4,4)) # Dummy dataset
    viewarr(data,index=1)

    EXAMPLE 3:
    from viewarr import *
    data=np.arange(64).reshape((4,4,4)) # Dummy dataset
    viewarr(data,index=1,idxnames=['ix','iy','iz'])

    EXAMPLE 4:
    from viewarr import *
    data=np.arange(64).reshape((4,4,4)) # Dummy dataset
    viewarr(data,index=1,idxnames=['x','y','z'],idxvals=[['a','b','c','d'],[-3,-1,1,3],[1.0,2.0,3.0,4.0]])

    EXAMPLE 5:
    from viewarr import *
    data1=np.arange(64).reshape((4,4,4)) # Dummy dataset
    data2=64-data1
    viewarr([data1,data2],index=1,idxnames=['x','y','z'],idxvals=[['a','b','c','d'],[-3,-1,1,3],[1.0,2.0,3.0,4.0]],ylabel=['Bla','adfsd'])
    """
    import matplotlib.pyplot as plt
    if type(data)==list:
        shape =  data[0].shape
        ndim  = len(shape)
    else:
        shape = data.shape
        ndim  = len(shape)
    assert index<ndim, "Index out of range"
    idxorder  = list(range(ndim))
    idxorder.pop(index)
    idxorder.append(index)
    if type(data)==list:
        datatrans = []
        for d in data:
            datatrans.append(d.transpose(idxorder))
        shapetrans = datatrans[0].shape
    else:
        datatrans  = data.transpose(idxorder)
        shapetrans = datatrans.shape
    def func(x,param,fixedpar={"datatrans":datatrans}):
        datatrans = fixedpar["datatrans"]
        if type(datatrans)==list:
            answer = []
            for dslice in datatrans:
                for i in range(len(param)):
                    dslice = dslice[param[i]]
                answer.append(dslice)
        else:
            dslice = datatrans
            for i in range(len(param)):
                dslice = dslice[param[i]]
            answer = dslice
        answer = np.array(answer)
        return answer
    params=[]
    for i in range(ndim-1):
        params.append(np.arange(shapetrans[i]))
    if x is None:
        if idxvals is None:
            x = np.arange(shapetrans[-1])
        else:
            x = np.array(idxvals[index])
    if ymin is None:
        if type(data)==list:
            ymin = []
            for d in data:
                ymin.append(d.min())
            ymin = np.array(ymin).min()
        else:
            ymin = data.min()
    if ymax is None:
        if type(data)==list:
            ymax = []
            for d in data:
                ymax.append(d.max())
            ymax = np.array(ymax).max()
        else:
            ymax = data.max()
    if idxvals is not None:
        paramsalt = []
        for i in range(ndim-1):
            paramsalt.append(idxvals[idxorder[i]])
    else:
        paramsalt = None
    fig = None
    ax  = None
    if idxnames is None:
        parnames = []
        for i in range(ndim-1):
            s = 'Parameter {}'.format(idxorder[i])
            parnames.append(s)
        xname    = 'Parameter {}'.format(index)
    else:
        parnames = []
        for i in range(ndim-1):
            parnames.append(idxnames[idxorder[i]]+" =")
        xname    = idxnames[index]
    fig = plt.figure()
    ax  = plt.axes(xlim=(x.min(),x.max()),ylim=(ymin,ymax))
    ax.set_xlabel(xname)
    if ylabel is not None:
        if type(ylabel)==list:
            label = r''
            glue  = ''
            for l in ylabel:
                label += glue+l
                glue = ', '
            ax.set_ylabel(label)
        else:
            ax.set_ylabel(ylabel)
    if type(data)==list:
        axmodel = []
        if ylabel is None:
            for i in range(len(datatrans)):
                axm0,  = ax.plot(x,x,label='{}'.format(i))
                axmodel.append(axm0)
        else:
            for i in range(len(datatrans)):
                axm0,  = ax.plot(x,x,label=ylabel[i])
                axmodel.append(axm0)
        ax.legend()
    else:
        axmodel = None
    interactive_plot(x, func, params, ymin=ymin, ymax=ymax, parnames=parnames, parunits=None, fig=fig, ax=ax, axmodel=axmodel, parstart=None, iparstart=None, plotbutton=False, fixedpar=None, returnipar=False, block=False, paramsalt=paramsalt, altformat=idxformat)

#
# Interactive plotting tool
#
# Copyright (c) 2018 C.P. Dullemond
# Free software under the standard MIT License
#
import numpy as np

def interactive_plot(x, func, params, ymin=None, ymax=None, parnames=None, parunits=None, \
                     fig=None, ax=None, axmodel=None, parstart=None, iparstart=None,      \
                     plotbutton=False, fixedpar=None, returnipar=False, block=False,      \
                     paramsalt=None, altformat='', img_x=None, img_y=None, img_func=None, \
                     img_im=None, **kwargs):
    """
    Plot the function func(x) with parameters given by the params
    list of lists. 

    ARGUMENTS:
      x          Array of x values
      func       Function func(x,params)
      params     List of parameters, but with each parameter value
                 here given as a list of possible values.

    OPTIONAL ARGUMENTS:
      ymin       Set vertical axis lower limit
      ymax       Set vertical axis upper limit
      parnames   Names of the params, e.g. ['A', 'omega']
                 If the parnames have an '=' sign (e.g. ['A = ', 'omega = '])
                 then the value of the parameters are written out.
      parunits   If set, a list of values by which the parameter values are divided
                 before being printed on the widget (only if parnames have '=').
                 It only affects the printing next to the sliders, and has no 
                 other effect.
      fig        A pre-existing figure
      ax         A pre-existing axis
      axmodel    If set, this is the plot style of the model
      parstart   If set, set the sliders initially close to these values
      iparstart  If set, set the slider index values initially to these values
                 (note: iparstart is an alternative to parstart)
      paramsalt  If set, then instead of the params values, the paramsalt values 
                 will be written after '=' (only if parnames is set, see above).
      returnipar If True, then return ipar
      block      If True, then wait until window is closed

    EXAMPLE 1 (Simplest example):
    from interactive_plot import *
    def func(x,param): return param[0]*np.sin(param[1]*x)
    x      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)] # Choices of parameter values
    interactive_plot(x, func, params, ymax=1., ymin=-1., parnames=['A = ','omega = '])

    EXAMPLE 1-a (With plotting button instead of automatic replot; useful for heavier models):
    from interactive_plot import *
    def func(x,param): return param[0]*np.sin(param[1]*x)
    x      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)] # Choices of parameter values
    interactive_plot(x, func, params, ymax=1., ymin=-1., parnames=['A = ','omega = '],plotbutton=True)

    EXAMPLE 1-b (Plotting the content of a pre-calculated 2-D array)
    from interactive_plot import *
    x       = np.linspace(0,2*np.pi,100)
    y_array = np.zeros((30,100))
    omega   = np.linspace(1,3.,30)
    for i in range(30): y_array[i,:] = np.sin(omega[i]*x)
    def func(x,param): return y_array[param[0],:]
    params  = [np.arange(30)] # Choices of parameter values
    interactive_plot(x, func, params)

    EXAMPLE 2 (Model fitting to data):
    import numpy as np
    import matplotlib.pyplot as plt
    from interactive_plot import *
    def func(x,param): return param[0]*np.sin(param[1]*x)
    x        = np.linspace(0,2*np.pi,100)
    data     = 0.5*np.sin(2.*x)*(1.0+0.6*np.random.normal(size=len(x)))
    fig      = plt.figure(1)
    ax       = plt.axes(xlim=(x.min(),x.max()),ylim=(-1.2,1.2))
    axd,     = ax.plot(x,data,'o',label='data')
    plt.xlabel('x [cm]')
    plt.ylabel('f [erg/s]')
    params   = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)] # Choices of parameter values
    parstart = [0.6,2.0]  # Initial guesses for parameters
    interactive_plot(x, func, params, parnames=['A = ','omega = '], fig=fig, ax=ax, label='model',parstart=parstart)
    ax.legend()
    plt.show()

    EXAMPLE 2-a (Model overplotting over an image):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from interactive_plot import *
    def func(x,param): return param[0]*np.sin(param[1]*x)
    x        = np.linspace(0,2*np.pi,100)
    image    = np.random.normal(size=(70,70)) # Make some image
    fig      = plt.figure(1)
    extent   = [x.min(),x.max(),-1.2,1.2]
    axd      = plt.imshow(image,extent=extent,cmap=cm.hot)
    ax       = plt.gca()
    plt.axis(extent)
    plt.xlabel('x [cm]')
    plt.ylabel('f [erg/s]')
    params   = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)] # Choices of parameter values
    parstart = [0.6,2.0]  # Initial guesses for parameters
    interactive_plot(x, func, params, parnames=['A = ','omega = '], fig=fig, ax=ax, label='model',parstart=parstart)
    ax.legend()
    plt.show()

    EXAMPLE 3 (Fitting two models simultaneously to data):
    import numpy as np
    import matplotlib.pyplot as plt
    from interactive_plot import *
    def func(x,param): return np.vstack((param[0]*np.sin(param[1]*x),param[0]*np.cos(param[1]*x)))
    x      = np.linspace(0,2*np.pi,100)
    data   = 0.5*np.sin(2.*x)*(1.0+0.6*np.random.normal(size=len(x)))
    fig    = plt.figure(1)
    ax     = plt.axes(xlim=(x.min(),x.max()),ylim=(-1.2,1.2))
    axd,   = ax.plot(x,data,'o',label='data')
    axm0,  = ax.plot(x,data,'--',label='sin')
    axm1,  = ax.plot(x,data,':',label='cos')
    axmodel= [axm0,axm1]
    plt.xlabel('x [cm]')
    plt.ylabel('f [erg/s]')
    params = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)]
    interactive_plot(x, func, params, parnames=['A = ','omega = '], fig=fig, ax=ax, axmodel=axmodel)
    ax.legend()
    plt.show()

    EXAMPLE 3-a (Fitting two models in two separate plots simultaneously):
    import numpy as np
    import matplotlib.pyplot as plt
    from interactive_plot import *
    def func(x,param): return np.vstack((param[0]*np.sin(param[1]*x),param[0]*np.cos(param[1]*x)))
    x         = np.linspace(0,2*np.pi,100)
    data      = 0.5*np.sin(2.*x)*(1.0+0.6*np.random.normal(size=len(x)))
    extent    = [x.min(),x.max(),-1.2,1.2]
    fig, axes = plt.subplots(ncols=2)
    axes[0].axis(extent)
    axes[1].axis(extent)
    axd0,  = axes[0].plot(x,data,'o',label='data')
    axm0,  = axes[0].plot(x,data,'--',label='sin')
    axd1,  = axes[1].plot(x,data,'o',label='data')
    axm1,  = axes[1].plot(x,data,':',label='cos')
    axmodel= [axm0,axm1]
    params = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)]
    interactive_plot(x, func, params, parnames=['A = ','omega = '], fig=fig, ax=0, axmodel=axmodel)
    plt.show()

    EXAMPLE 4: (passing additional fixed parameters to function):
    from interactive_plot import *
    def func(x,param,fixedpar={}): return param[0]*np.sin(param[1]*x)+fixedpar['offset']
    x      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(1.,3.,30)] # Choices of parameter values
    interactive_plot(x, func, params, ymax=1., ymin=-1., parnames=['A = ','omega = '],fixedpar={'offset':0.6})

    EXAMPLE 5: (Interactive image, i.e. 2D slice from a higher-dimensional data box)
    import numpy as np
    from interactive_plot import *
    from matplotlib import cm
    from matplotlib import colors
    import matplotlib.pyplot as plt
    from matplotlib.image import NonUniformImage
    x        = np.linspace(-1,1,20)
    y        = np.linspace(-1,1,30)
    z        = np.linspace(0,1,25)
    xx,yy,zz = np.meshgrid(x,y,z,indexing='ij')
    rr       = np.sqrt(xx**2+yy**2)
    f        = np.sin(xx*2*np.pi)*yy*(1-zz)+np.cos(2*np.pi*rr)*zz
    norm     = colors.Normalize(vmin=f.min(),vmax=f.max())
    cmap     = cm.hot
    fig,ax   = plt.subplots()
    im       = NonUniformImage(ax,interpolation='nearest',cmap=cmap,norm=norm)
    im.set_data(x,y,f[:,:,0].T)
    ax.images.append(im)
    ax.set_xlim((x[0]-0.5*(x[1]-x[0]),x[-1]+0.5*(x[-1]-x[-2])))
    ax.set_ylim((y[0]-0.5*(y[1]-y[0]),y[-1]+0.5*(y[-1]-y[-2])))
    cbar=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax)
    cbar.set_label(r'$T\;[\mathrm{K}]$')
    def img_func(param,fixedpar={}): return fixedpar['f'][:,:,param[0]]
    params = [np.arange(25)] # Choices of parameter values
    fixedpar = {}
    fixedpar["f"]=f
    interactive_plot(None, None, params, fixedpar=fixedpar,       \
                     img_x=x,img_y=y,img_func=img_func,img_im=im, \
                     fig=fig,ax=ax)
    
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons

    # Compute spacing of plot, sliders and button
    hslider  = 0.03
    nslidrscl= 6
    if(len(params)>nslidrscl):
        hslider *= float(nslidrscl)/len(params)
    dyslider = hslider*(4./3.)
    xslider  = 0.3
    wslider  = 0.3
    hbutton  = 0.06
    wbutton  = 0.15
    xbutton  = 0.3
    dybutton = hbutton+0.01
    panelbot = 0.0
    controlh = panelbot + len(params)*dyslider
    if plotbutton: controlh += dybutton
    controltop = panelbot + controlh
    bmargin  = 0.15
    
    # generate figure
    if fig is None: fig = plt.figure()
    fig.subplots_adjust(top=0.95,bottom=controltop+bmargin)

    # Set the initial values
    indexinit = np.zeros(len(params),dtype=int)
    if parstart is not None:
        for i in range(len(params)):
            if parstart[i] in params[i]:
                idx = np.where(np.array(params[i])==parstart[i])[0]
                if len(idx)>0:
                    indexinit[i] = idx[0]
            else:
                if params[i][-1]>params[i][0]:
                    idx = np.where(np.array(params[i])<parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[-1]
                else:
                    idx = np.where(np.array(params[i])>parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[0]
    if iparstart is not None:
        indexinit[:] = iparstart[:]

    # select first image
    par = []
    for i in range(len(params)):
        par.append(params[i][indexinit[i]])
    xmin = None
    xmax = None
    if x is not None:
        if xmin is None: xmin = x.min()
        if xmax is None: xmax = x.max()
    if func is not None:
        if fixedpar is not None:
            f = func(x,par,fixedpar=fixedpar)
        else:
            f = func(x,par)
        # set range
        if ymin is None: ymin = f.min()
        if ymax is None: ymax = f.max()
    if img_y is not None:
        if ymin is None: ymin = img_y.min()
        if ymax is None: ymax = img_y.max()
        if ymin>img_y.min(): ymin = img_y.min()
        if ymax<img_y.max(): ymax = img_y.max()
    if img_x is not None:
        if xmin is None: xmin = img_x.min()
        if xmax is None: xmax = img_x.max()
        if xmin>img_x.min(): xmin = img_x.min()
        if xmax<img_x.max(): xmax = img_x.max()
    assert (xmin is not None) or (xmax is not None) , 'Error: x undefined'
    assert (ymin is not None) or (ymay is not None) , 'Error: y undefined'
    
    # display function(s)
    if ax is None:      ax       = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
    if axmodel is None:
        if func is not None:
            if len(f.shape)==1:
                # Normal case: a single model function
                axmodel, = ax.plot(x,f,**kwargs)
            else:
                # Special case: multiple model functions: f[imodel,:]
                assert len(f.shape)==2, 'Model returns array with more than 2 dimensions. No idea what to do.'
                axmodel = []
                for i in range(f.shape[0]):
                    axm, = ax.plot(x,f[i,:],**kwargs)
                    axmodel.append(axm)
            
    sliders = []
    for i in range(len(params)):
    
        # define slider
        axcolor = 'lightgoldenrodyellow'
        axs = fig.add_axes([xslider, controltop-i*dyslider, xslider+wslider, hslider], facecolor=axcolor)

        if parnames is not None:
            name = parnames[i]
        else:
            name = 'Parameter {0:d}'.format(i)

        slider = Slider(axs, name, 0, len(params[i]) - 1,
                    valinit=indexinit[i], valfmt='%i')
        sliders.append(slider)

    if plotbutton:
        axb = fig.add_axes([xbutton, panelbot+0.2*hbutton, xbutton+wbutton, hbutton])
        pbutton = Button(axb,'Plot')
    else:
        pbutton = None

    class callbackplot(object):
        def __init__(self,x,func,params,sliders,pbutton=None,fixedpar=None,ipar=None, \
                     img_x=None, img_y=None, img_func=None, img_im=None):
            self.x        = x
            self.func     = func
            self.params   = params
            self.sliders  = sliders
            self.pbutton  = pbutton
            self.fixedpar = fixedpar
            self.parunits = parunits
            self.paramsalt= paramsalt
            self.altformat= altformat
            self.img_x    = img_x
            self.img_y    = img_y
            self.img_func = img_func
            self.img_im   = img_im
            self.closed   = False
            if ipar is None:
                self.ipar = np.zeros(len(sliders),dtype=int)
            else:
                self.ipar = ipar
        def handle_close(self,event):
            self.closed   = True
        def myreadsliders(self):
            for isl in range(len(self.sliders)):
                ind = int(self.sliders[isl].val)
                self.ipar[isl]=ind
            par = []
            for i in range(len(self.ipar)):
                ip = self.ipar[i]
                value = self.params[i][ip]
                par.append(value)
                name = self.sliders[i].label.get_text()
                if '=' in name:
                    namebase = name.split('=')[0]
                    if self.paramsalt is not None:
                        vls  = "{0:" + self.altformat + "}"
                        name = namebase + "= " + vls.format(self.paramsalt[i][ip])
                    else:
                        if self.parunits is not None:
                            valunit = self.parunits[i]
                        else:
                            valunit = 1.0
                        name = namebase + "= {0:13.6e}".format(value/valunit)
                    self.sliders[i].label.set_text(name)
            return par
        def myreplot(self,par):
            if self.x is not None and self.func is not None:
                x = self.x
                if self.fixedpar is not None:
                    f = self.func(x,par,fixedpar=self.fixedpar)
                else:
                    f = self.func(x,par)
                if len(f.shape)==1:
                    axmodel.set_data(x,f)
                else:
                    for i in range(f.shape[0]):
                        axmodel[i].set_data(x,f[i,:])
            if self.img_x is not None and self.img_y is not None and self.img_func is not None and self.img_im is not None:
                x = self.img_x
                y = self.img_y
                if self.fixedpar is not None:
                    z = self.img_func(par,fixedpar=self.fixedpar)
                else:
                    z = self.img_func(par)
                self.img_im.set_data(x,y,z.T)
            plt.draw()
        def mysupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is None: self.myreplot(par)
        def mybupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is not None: self.pbutton.label.set_text('Computing...')
            plt.pause(0.01)
            self.myreplot(par)
            if self.pbutton is not None: self.pbutton.label.set_text('Plot')

    mcb = callbackplot(x,func,params,sliders,pbutton=pbutton,fixedpar=fixedpar,ipar=indexinit, \
                       img_x=img_x,img_y=img_y,img_func=img_func,img_im=img_im)

    mcb.mybupdate(0)

    if plotbutton:
        pbutton.on_clicked(mcb.mybupdate)
    for s in sliders:
        s.on_changed(mcb.mysupdate)

    fig._mycallback    = mcb

    if block:
        plt.show(block=True)
    if returnipar:
        return mcb.ipar
        

def interactive_curve(t, func, params, xmin=None, xmax=None, ymin=None, ymax=None, parnames=None, parunits=None, fig=None, ax=None, axmodel=None, parstart=None, iparstart=None, plotbutton=False, fixedpar=None, returnipar=False, block=False, **kwargs):
    """
    Plot the 2-D curve x,y = func(t) with parameters given by the params
    list of lists. 

    ARGUMENTS:
      t          Array of t values
      func       Function func(x,params)
      params     List of parameters, but with each parameter value
                 here given as a list of possible values.

    OPTIONAL ARGUMENTS:
      xmin       Set horizontal axis lower limit
      xmax       Set horizontal axis upper limit
      ymin       Set vertical axis lower limit
      ymax       Set vertical axis upper limit
      parnames   Names of the params, e.g. ['A', 'omega']
                 If the parnames have an '=' sign (e.g. ['A = ', 'omega = '])
                 then the value of the parameters are written out.
      parunits   If set, a list of values by which the parameter values are divided
                 before being printed on the widget (only if parnames have '=').
                 It only affects the printing next to the sliders, and has no 
                 other effect.
      fig        A pre-existing figure
      ax         A pre-existing axis
      parstart   If set, set the sliders initially close to these values
      iparstart  If set, set the slider index values initially to these values
                 (note: iparstart is an alternative to parstart)
      returnipar If True, then return ipar
      block      If True, then wait until window is closed

    EXAMPLE 1 (one ellipse):
    from interactive_plot import *
    def func(t,param): 
        x = param[0]*np.cos(t)
        y = param[1]*np.sin(t)
        csw = np.cos(param[2])
        snw = np.sin(param[2])
        return csw*x-snw*y,snw*x+csw*y
    t      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(0.1,1.,30),np.linspace(0.,np.pi,30)]
    interactive_curve(t, func, params, xmax=1., xmin=-1., ymax=1., ymin=-1., parnames=['Ax = ','Ay = ','omega = '],iparstart=[10,15,12])

    EXAMPLE 1-a (With plotting button instead of automatic replot; useful for heavier models):
    from interactive_plot import *
    def func(t,param): 
        x = param[0]*np.cos(t)
        y = param[1]*np.sin(t)
        csw = np.cos(param[2])
        snw = np.sin(param[2])
        return csw*x-snw*y,snw*x+csw*y
    t      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(0.1,1.,30),np.linspace(0.,np.pi,30)]
    interactive_curve(t, func, params, xmax=1., xmin=-1., ymax=1., ymin=-1., parnames=['Ax = ','Ay = ','omega = '],iparstart=[10,15,12],plotbutton=True)

    EXAMPLE 2 (two ellipses):
    import numpy as np
    import matplotlib.pyplot as plt
    from interactive_plot import *
    def func(t,param): 
        x = param[0]*np.cos(t)
        y = param[1]*np.sin(t)
        csw = np.cos(param[2])
        snw = np.sin(param[2])
        return np.vstack((csw*x-snw*y,-csw*x-snw*y)),np.vstack((snw*x+csw*y,snw*x+csw*y))
    t      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(0.1,1.,30),np.linspace(0.,np.pi,30)]
    fig    = plt.figure(1)
    ax     = plt.axes(xlim=(-1.2,1.2),ylim=(-1.2,1.2))
    x,y    = func(t,[1.,1.,1.])
    axm0,  = ax.plot(x[0,:],y[0,:],'--',label='left')
    axm1,  = ax.plot(x[1,:],y[1,:],':',label='right')
    axmodel= [axm0,axm1]
    interactive_curve(t, func, params, xmax=1., xmin=-1., ymax=1., ymin=-1., parnames=['Ax = ','Ay = ','omega = '],iparstart=[10,15,12], fig=fig, ax=ax, axmodel=axmodel)

    EXAMPLE 3 (as example 2, but now each ellipse in its own panel):
    import numpy as np
    import matplotlib.pyplot as plt
    from interactive_plot import *
    def func(t,param): 
        x = param[0]*np.cos(t)
        y = param[1]*np.sin(t)
        csw = np.cos(param[2])
        snw = np.sin(param[2])
        return np.vstack((csw*x-snw*y,-csw*x-snw*y)),np.vstack((snw*x+csw*y,snw*x+csw*y))
    t      = np.linspace(0,2*np.pi,100)
    params = [np.linspace(0.1,1.,30),np.linspace(0.1,1.,30),np.linspace(0.,np.pi,30)]
    fig, axes = plt.subplots(nrows=2)
    axes[0].set_xlim((-1.2,1.2))
    axes[0].set_ylim((-1.2,1.2))
    axes[1].set_xlim((-1.2,1.2))
    axes[1].set_ylim((-0.8,0.8))
    x,y    = func(t,[1.,1.,1.])
    axm0,  = axes[0].plot(x[0,:],y[0,:],'--',label='left')
    axm1,  = axes[1].plot(x[1,:],y[1,:],':',label='right')
    axmodel= [axm0,axm1]
    interactive_curve(t, func, params, xmax=1., xmin=-1., ymax=1., ymin=-1., parnames=['Ax = ','Ay = ','omega = '],iparstart=[10,15,12], fig=fig, ax=axes[0], axmodel=axmodel)
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons

    # Compute spacing of plot, sliders and button
    hslider  = 0.03
    nslidrscl= 6
    if(len(params)>nslidrscl):
        hslider *= float(nslidrscl)/len(params)
    dyslider = hslider*(4./3.)
    xslider  = 0.3
    wslider  = 0.3
    hbutton  = 0.06
    wbutton  = 0.15
    xbutton  = 0.3
    dybutton = hbutton+0.01
    panelbot = 0.0
    controlh = panelbot + len(params)*dyslider
    if plotbutton: controlh += dybutton
    controltop = panelbot + controlh
    bmargin  = 0.15
    
    # generate figure
    if fig is None: fig = plt.figure()
    fig.subplots_adjust(top=0.95,bottom=controltop+bmargin)

    # Set the initial values
    indexinit = np.zeros(len(params),dtype=int)
    if parstart is not None:
        for i in range(len(params)):
            if parstart[i] in params[i]:
                idx = np.where(np.array(params[i])==parstart[i])[0]
                if len(idx)>0:
                    indexinit[i] = idx[0]
            else:
                if params[i][-1]>params[i][0]:
                    idx = np.where(np.array(params[i])<parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[-1]
                else:
                    idx = np.where(np.array(params[i])>parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[0]
    if iparstart is not None:
        indexinit[:] = iparstart[:]

    # select first image
    par = []
    for i in range(len(params)):
        par.append(params[i][indexinit[i]])
    if fixedpar is not None:
        x, y = func(t,par,fixedpar=fixedpar)
    else:
        x, y = func(t,par)

    # set range
    if xmin is None: xmin = x.min()
    if xmax is None: xmax = x.max()
    if ymin is None: ymin = y.min()
    if ymax is None: ymax = y.max()
    
    # display function
    if ax is None: ax   = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
    if axmodel is None:
        if len(x.shape)==1:
            # Normal case: a single model function
            assert len(x.shape)==1, 'Cannot have multiple y and single x'
            axmodel, = ax.plot(x,y,**kwargs)
        else:
            # Special case: multiple model functions: f[imodel,:]
            assert len(x.shape)==2, 'Model returns array with more than 2 dimensions. No idea what to do.'
            assert len(y.shape)==2, 'Cannot have multiple x and single y'
            axmodel = []
            for i in range(x.shape[0]):
                axm, = ax.plot(x[i,:],y[i,:],**kwargs)
                axmodel.append(axm)
    
    sliders = []
    for i in range(len(params)):
    
        # define slider
        axcolor = 'lightgoldenrodyellow'
        axs = fig.add_axes([xslider, controltop-i*dyslider, xslider+wslider, hslider], facecolor=axcolor)

        if parnames is not None:
            name = parnames[i]
        else:
            name = 'Parameter {0:d}'.format(i)
            
        slider = Slider(axs, name, 0, len(params[i]) - 1,
                    valinit=indexinit[i], valfmt='%i')
        sliders.append(slider)

    if plotbutton:
        axb = fig.add_axes([xbutton, panelbot+0.2*hbutton, xbutton+wbutton, hbutton])
        pbutton = Button(axb,'Plot')
    else:
        pbutton = None

    class callbackcurve(object):
        def __init__(self,t,func,params,sliders,pbutton=None,fixedpar=None,ipar=None):
            self.t        = t
            self.func     = func
            self.params   = params
            self.sliders  = sliders
            self.pbutton  = pbutton
            self.fixedpar = fixedpar
            self.parunits = parunits
            self.closed   = False
            if ipar is None:
                self.ipar = np.zeros(len(sliders),dtype=int)
            else:
                self.ipar = ipar
        def handle_close(self,event):
            self.closed   = True
        def myreadsliders(self):
            for isl in range(len(self.sliders)):
                ind = int(self.sliders[isl].val)
                self.ipar[isl]=ind
            par = []
            for i in range(len(self.ipar)):
                ip = self.ipar[i]
                value = self.params[i][ip]
                par.append(value)
                name = self.sliders[i].label.get_text()
                if '=' in name:
                    namebase = name.split('=')[0]
                    if self.parunits is not None:
                        valunit = self.parunits[i]
                    else:
                        valunit = 1.0
                    name = namebase + "= {0:13.6e}".format(value/valunit)
                    self.sliders[i].label.set_text(name)
            return par
        def myreplot(self,par):
            t = self.t
            if self.fixedpar is not None:
                x,y = self.func(t,par,fixedpar=self.fixedpar)
            else:
                x,y = self.func(t,par)
            if len(x.shape)==1:
                axmodel.set_data(x,y)
            else:
                for i in range(x.shape[0]):
                    axmodel[i].set_data(x[i,:],y[i,:])
            plt.draw()
        def mysupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is None: self.myreplot(par)
        def mybupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is not None: self.pbutton.label.set_text('Computing...')
            plt.pause(0.01)
            self.myreplot(par)
            if self.pbutton is not None: self.pbutton.label.set_text('Plot')

    mcb = callbackcurve(t,func,params,sliders,pbutton=pbutton,fixedpar=fixedpar,ipar=indexinit)
            
    mcb.mybupdate(0)
        
    if plotbutton:
        pbutton.on_clicked(mcb.mybupdate)
    for s in sliders:
        s.on_changed(mcb.mysupdate)

    fig._mycallback    = mcb
    
    if block:
        plt.show(block=True)
    if returnipar:
        return mcb.ipar
