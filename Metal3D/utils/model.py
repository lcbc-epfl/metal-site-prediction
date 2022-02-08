import warnings 

import torch
import torch.nn as nn
import torch.nn.functional as F

class Model(nn.Module):
    """ Model with same padding
    Conv5 uses a large filter size to aggregate the features from the whole box """
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Conv3d(8, 32, 3, padding='same') 
        self.conv2 = nn.Conv3d(32, 64, 4,padding='same')  
        self.conv3 = nn.Conv3d(64, 80, 3,padding='same')
        self.conv4 = nn.Conv3d(80, 20, 3,padding='same')
        self.conv5 = nn.Conv3d(20, 20, 16,padding='same')
        self.conv6 = nn.Conv3d(20, 16, 3,padding='same')
        self.conv7 = nn.Conv3d(16, 1, 3,padding='same')
        self.dropout1 = nn.Dropout(0.2) 
        self.dropout2 = nn.Dropout(0.2)
        
    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x= self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)

        x = self.conv4(x)
        x = F.relu(x)

        x = self.conv5(x)
        x = F.relu(x)
        x = self.dropout1(x)
        x = self.conv6(x)
        x = F.relu(x)

        x = self.conv7(x)
        x = torch.sigmoid(x)
        return x