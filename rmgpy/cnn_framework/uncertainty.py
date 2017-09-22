import keras.backend as K
import numpy as np
from keras.engine.topology import Layer 
from keras.models import Model
import time

class EnsembleModel(Model):
    def __init__(self, seeds=[], **kwargs):
        super(EnsembleModel, self).__init__(**kwargs)
        self.grad_model = None
        self.seeds = seeds
        self.weight_generators = None
        
    def gen_mask(self, seed):
        np.random.seed(seed)
        for layer in self.layers:
            if 'gen_mask' in dir(layer):
                layer.gen_mask()
        
    def set_weight_generators(self):
        self.weight_generators = [np.random.RandomState(seed) for seed in self.seeds]

    def train_on_batch(self, x, y, **kwargs):   
        n_data = len(y)     
        losses = []
        
        for j in range(len(self.seeds)):
            self.gen_mask(self.seeds[j])
            weight = self.weight_generators[j].dirichlet(np.ones(n_data))*n_data               
            loss = super(EnsembleModel,self).train_on_batch(x, y, sample_weight=weight) 
            losses.append(loss) 
        
        return np.mean(losses)

    def test_on_batch(self, x, y, **kwargs):   
        n_data = len(y)     
        losses = []
        
        for j in range(len(self.seeds)):
            self.gen_mask(self.seeds[j])
            weight = self.weight_generators[j].dirichlet(np.ones(n_data))*n_data              
            loss = super(EnsembleModel,self).test_on_batch(x, y, sample_weight=weight) 
            losses.append(loss) 
        
        return np.mean(losses)
    
    # def fit(self, x, y, n_epoch=25, batch_size=50):
    #     n_data = len(x)
    #     n_batch = int(np.ceil(float(n_data) / batch_size))
    #     seeds = self.seeds
    #     n_set = len(self.seeds)
    #     training_order = range(n_data)
    #     np.random.shuffle(training_order)  
    #     
    #     for i in range(n_epoch):
    #         for j in range(n_set):
    #             set_training_start = time.time()
    #             RandomMask.gen_masks(seeds[j])
    #             np.random.seed(seeds[j])
    #             weights = np.random.dirichlet(np.ones(n_data), 1)[0]            
    #             losses = []
    #             for k in range(n_batch):
    #                 
    #                 start = k*batch_size
    #                 end = min(start + batch_size, n_data)
    #                 
    #                 x_batch = x[training_order[start:end]]
    #                 y_batch = y[training_order[start:end]]
    #                 weights_batch = weights[training_order[start:end]]            
    #                                     
    #                 loss = self.train_on_batch(x_batch, y_batch, weights_batch) 
    #                 losses.append(loss)
    #             set_training_end = time.time()   
    #                           
    #             print 'epoch {} set {} loss {} time {}'.format(i, j, np.mean(losses), set_training_end - set_training_start)      

    def predict(self, x, **kwargs):
        Y = []
        for j in range(len(self.seeds)):
            self.gen_mask(self.seeds[j])
            Y += [super(EnsembleModel,self).predict(x, **kwargs)] 
        Y_avg = np.mean(Y,axis=0)
    #    Y_var = 1.0/tou + np.var(Y,axis=0)
        #Y_var = np.var(Y,axis=0)
        return Y_avg  
    
    def variance(self, x, **kwargs):
        Y = []
        for j in range(len(self.seeds)):
            self.gen_mask(self.seeds[j])
            Y += [super(EnsembleModel,self).predict(x, **kwargs)] 
    #    Y_avg = np.mean(Y,axis=0)
    #    Y_var = 1.0/tou + np.var(Y,axis=0)
        Y_var = np.var(Y,axis=0)
        return Y_var       
    
    # def variance_grad(self, x):
    #     if self.grad_model is None:
    #         # Currently, we assume only one input, one output
    #         grad = K.gradients(self.outputs[0], self.inputs[0])
    #         self.grad_model = K.function(self.inputs, self.outputs+grad)
    #     
    #     Y = []
    #     dY = []
    #     shape = (-1, self.inputs[0].get_shape()[1])
    #     x = x.reshape(shape)
    #     for j in range(len(self.seeds)):
    #         self.gen_mask(self.seeds[j])
    #         #results = self.grad_model([x.reshape((-1,1))])
    #         results = self.grad_model([x])
    #         Y.append(results[0])
    #         dY.append(results[1])
    #          
    #     Y_avg = np.mean(Y,axis=0)
    #     dY_avg = np.mean(dY,axis=0)
    #     dY_var = 0.0
    #     for i in range(len(seeds)):
    #         dY_var += 2.0*(Y[i]-Y_avg)*(dY[i]-dY_avg) 
    # 
    #     return dY_var
    
    def get_config(self):
        config = super(EnsembleModel, self).get_config()
        config['seeds'] = self.seeds
        return config
    
    @classmethod
    def from_config(cls, config, custom_objects=None):
        model = super(EnsembleModel, cls).from_config(config, custom_objects=custom_objects)
        model.seeds = config.get('seeds')
        return model
        
class RandomMask(Layer):
    """Applies Mask to the input.
    """
    # masks = []
    # rates = []
    # 
    # @classmethod
    # def gen_masks(cls, seed):
        # np.random.seed(seed)
        # for i, mask in enumerate(cls.masks):
        #     retain_prob = 1.0 - cls.rates[i]
        #     size = K.int_shape(mask)
        #     K.set_value(mask,np.random.binomial(n=1,p=retain_prob,size=size).astype(np.float32))
        # 
#    @classmethod
#    def initialize(cls, n_set, seed=None):
#        cls.n_set = n_set
#        if seed is not None:
#            np.random.seed(seed)
#        cls.seeds = np.random.randint(0, 10e8, n_set)
    
    def __init__(self, dropout_rate, **kwargs):
        self.dropout_rate = dropout_rate       
        super(RandomMask, self).__init__(**kwargs)

    def call(self, x, **kwargs):
        size = K.int_shape(x)[1:]
        self.mask = K.variable(np.ones(shape=size,dtype=np.float32))
        x *= self.mask
        return x
    
    def gen_mask(self):
        retain_prob = 1.0 - self.dropout_rate
        size = K.int_shape(self.mask)
        K.set_value(self.mask,np.random.binomial(n=1,p=retain_prob,size=size).astype(np.float32))
    
    def get_config(self):
        config = {'dropout_rate': self.dropout_rate}
        base_config = super(RandomMask, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))
