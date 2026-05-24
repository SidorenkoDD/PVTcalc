class Conditions:
    '''
    Класс для хранения термобарических условий
    '''
    def __init__(self, p_bar, t_c):
        self.p = p_bar
        self.t = t_c + 273.14

class StandardConditions:
    '''
    Class for standard conditions
    '''
    def __init__(self):
        self.p = 0.101325
        self.t = 293.14