from numpy import *
from numpy.random import uniform,randint
from scipy.special import factorial,binom
from matplotlib.pyplot import *


####################### function definitions ############################

def collision_probability_theory(n_people,n_days=365) :
    """
    returns analytical collision probability for the birthday problem.
    
    --- paramters ---
    n_people : <ndarray> type:int
        array of total number of people in the given room
    n_days : <ndarray> type:int
        array of total number of possible days to have a birthday on

    --- returns ---
    collision_probability : <ndarray> type:float
        probabilities of having the same birthday in the same room
    
    Mathis, Frank H. 1991. A Generalized Birthday Problem. SIAM Review.
    Society for Industrial and Applied Mathematics. 33 (2): 265-270
    """

    unique_probability = factorial(n_people,exact=True) * binom(n_days,n_people) / float(n_days)**n_people
    return 1-unique_probability


def collision_probability_experiment(n_people,n_days=365,n_rooms=1000) :
    """
    generates rooms of people from uniform distributions and returns
    calculated collision probability for the birthday problem.
    
    --- paramters ---
    n_people : <ndarray> type:int
        array of total number of people in the given room
    n_days : <int>
        array of total number of possible days to have a birthday on
    n_rooms : <int>
        number of rooms to calculate probabilities over
        
    --- returns ---
    collision_probability : <ndarray> type:float
        probabilities of having the same birthday in the same room
    
    """
    
    # create rooms of people
    rooms = randint(1,n_days+1,size=(n_rooms,n_people))
    
    # check if people are unique in each room
    is_unique = array([ len(unique(people)) is n_people for people in rooms ])
    
    # proability of not being unique
    return 1-mean(is_unique)


# upgrades functions to accept ndarrays, numpy-style ;P
collision_probability_theory = vectorize(collision_probability_theory)
collision_probability_experiment = vectorize(collision_probability_experiment)
    

####################### main program ############################
def main() :

    # calculations for theory and numerical investigations
    n_people = linspace(0,100,50).astype(int)
    theory = collision_probability_theory(n_people)
    experiment = collision_probability_experiment(n_people)

    # compare the two in a plot
    figure(figsize=(10,10))
    plot(n_people,theory,'k',label='theory')
    plot(n_people,experiment,'r.',label='experiment')

    xlabel(r'Number of people, $N$',fontsize=16)
    ylabel(r'Probability of Collision, $p(N)$',fontsize=16)
    legend(fontsize=16)

    show()
    
    
# execute main program if this file run directly
if __name__== "__main__":
    main()
