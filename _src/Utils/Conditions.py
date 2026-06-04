class Conditions:
    '''
    Класс для хранения термобарических условий
    '''
    def __init__(self, p, t):
        self.p = p
        self.t = t + 273.14

class StandardConditions:
    '''
    Class for standard conditions
    '''
    def __init__(self):
        self.p = 1.01325
        self.t = 293.14