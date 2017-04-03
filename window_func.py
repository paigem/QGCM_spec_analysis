import numpy as np

# Function for 2d spatial window
def constructoverlappingwindows_2d_TUKEY(len1,len2):
	
	window1 = constructoverlappingwindows_1d_TUKEY(len1)
	window2 = constructoverlappingwindows_1d_TUKEY(len2)
	
	print 'window1.shape=',window1.shape,np.transpose(window1).shape
	
	# Construct the overlapping windows
	window = np.tile(window2,(len1,1))*np.transpose(np.tile(window1,(len2,1)))
	
	del window1,window2
	return window
	
# Function to construct 1d temporal window	
def constructoverlappingwindows_1d_TUKEY(len):

	# len should be multiple of 10 <--!!!!!CHECK WHY!!!!!!
	window = tukey(len,alpha=0.2,sym=True)
	
	return window
	
	
	
def tukey(M,alpha,sym):
    # r"""Return a Tukey window, also known as a tapered cosine window.
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')

    if alpha <= 0:
        return np.ones(M, 'd')
    elif alpha >= 1.0:
        return hann(M, sym=sym)

    odd = M % 2
    if not sym and not odd:
        M = M + 1

    n = np.arange(0, M)
    width = int(np.floor(alpha*(M-1)/2.0))
    n1 = n[0:width+1]
    n2 = n[width+1:M-width-1]
    n3 = n[M-width-1:]

    w1 = 0.5 * (1 + np.cos(np.pi * (-1 + 2.0*n1/alpha/(M-1))))
    w2 = np.ones(n2.shape)
    w3 = 0.5 * (1 + np.cos(np.pi * (-2.0/alpha + 1 + 2.0*n3/alpha/(M-1))))

    w = np.concatenate((w1, w2, w3))

    if not sym and not odd:
        w = w[:-1]
    del w1,w2,w3,n,M
    return w

def main(var,spacetime):

	# 2d spatial window
	window_2d = constructoverlappingwindows_2d_TUKEY(var.shape[0],var.shape[1])
		
	# 1d spatial window
	window_1d = constructoverlappingwindows_1d_TUKEY(var.shape[2])

	if spacetime == 'spacetime':
		print 'Making both time and space window'

		# Multiply var by windows in the relevant dimensions
		var = var*np.rollaxis(np.tile(window_2d,(var.shape[2],1,1)),0,3)*np.tile(window_1d,(var.shape[0],var.shape[1],1))
	
		print 'var.shape=',var.shape
		
	elif spacetime == 'space':
		print 'Making window in space'

		# Multiply var by windows in the relevant dimensions
		var = var*np.rollaxis(np.tile(window_2d,(var.shape[2],1,1)),0,3)

	elif spacetime == 'time':
		print 'Making window in time'

		# Multiply var by windows in the relevant dimensions
		var = var*np.tile(window_1d,(var.shape[0],var.shape[1],1))


	return var





