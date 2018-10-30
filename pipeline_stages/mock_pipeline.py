'''
Code to make a simple run of the pipeline
'''

class MockPipe():
    '''
    class to make a simple of run all stages of the pipeline
    '''
    def __init__(self):
        pass

    def stage1(self):
        '''
        First stage of the CLMM pipeline where we create the galaxy cluster
        objects!
    
        The inputs are galaxy cluster data, and the outputs are the
        populated galaxy cluster objects.
    
        The pipeline loops through input galaxy cluster data, collecting
        the information into galaxy cluster objects.
    
        It can do this in parallel if needed.  We might want to move some
        of the functionality here (e.g. the I/O) into a general parent
        class.
    
        '''
        print('*** stage 1 ***')

    def stage2(self):
        '''
        Second stage of the CLMM pipeline where we use the inferrer to sample
        likelihood of parameters on each cluster object.
    
        The inputs are the populated galaxy cluster objects, the outputs
        are the chains.
    
        The pipeline loops through input galaxy cluster objects, using the
        inference manager to collect relevant data for the Inferrer object. 
    
        It can do this in parallel if needed.  We might want to move some
        of the functionality here (e.g. the I/O) into a general parent
        class.
    
        '''
        print('*** stage 2 ***')

    def stage3(self):
        '''
        Third stage of the pipeline where we define the binning of for
        importance sampling.
    
        The inputs are galaxy cluster data with individual chains per
        cluster, and the outputs are the collections of cluster chains.
    
        The pipeline loops through galaxy cluster objects and identifies
        which collection they belong to.  Collections can be defined in
        richness bins, true mass bins, etc.
    
        It can do this in parallel if needed.  We might want to move some
        of the functionality here (e.g. the I/O) into a general parent
        class.
    
        '''
        print('*** stage 3 ***')

    def stage4(self):
        '''
        Fourth stage of the pipeline where we importance sample chains in
        each collection.
    
        The inputs are collections of chains, and the outputs are the
        results of importance sampling.
    
        The pipeline loops through the collections, and performs
        importance sampling.  Collections were defined in richness bins,
        true mass bins, etc. in the previous pipeline stage.
    
        It can do this in parallel if needed.  We might want to move some
        of the functionality here (e.g. the I/O) into a general parent
        class.
    
        '''
        print('*** stage 4 ***')

    def stage5(self):
        '''
        Fifth stage of the pipeline where we summarize the results of
        importance sampled chains per cluster collection.
    
        The inputs are importance sampled outputs per collection, and the
        outputs are the summary statistic(s). e.g. most likely mass with
        errorbars, concentration, noise levels, etc.
    
        The pipeline loops through the collections' outputs, and
        summarizes.
    
        It can do this in parallel if needed.  We might want to move some
        of the functionality here (e.g. the I/O) into a general parent
        class.
    
        '''
        print('*** stage 5 ***')

if __name__=='__main__':

    mp = MockPipe()

    mp.stage1()
