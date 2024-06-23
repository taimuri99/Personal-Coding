from numpy import *
from numpy.random import uniform,randint,choice,normal
from scipy.special import factorial,binom
from string import ascii_uppercase
from matplotlib.pyplot import *


####################### class definitions ############################
class Student(object) :
    """
    Object that stores all you need to know about a set of students.
    The student object behaves like a little bit like a vector.
    
    Adding two students gives another student whos height and birthday
    is the sum. Initials are added to the set of initials to keep track
    of how many students were operated on together.
    
    Adding a student to a number adds that number to height and birthday
    while the initials are keept the same. Dividing a student by a
    number divides their height and birthday by the same number and
    initials still the same.
    
    --- attributes ---
    
    self.initials : set(...<str>...)
        set of initials of the student(s)
        
    self.height : <float>
        their height in cm
        
    self.birthday : <int>
        month/day they were born

    """
    def __init__(self,initials,height,birthday) :
        """this is the constructor, called upon instantiation"""

        # forces correct types for each attribute
        self.initials = set(str(initials)) if type(initials) is not set else initials
        self.height = float(height)
        self.birthday = int(birthday)
        
    def __add__(self,instance) :
        """this is the addition operator, called upon addition of two students"""
        
        if  type(instance) is Student : # rules for adding two students
            
            initials = self.initials.union(instance.initials) 
            height = self.height + instance.height 
            birthday = self.birthday + instance.birthday
            
        elif type(instance) is float or type(instance) is int : # rules for adding a student to number
            
            initials = self.initials
            height = self.height + instance
            birthday = self.birthday + instance
            
        else : # raise error for anything else
            raise Exception('addition not supported for ',type(instance))

        return Student(initials,height,birthday)

    
    def __div__(self,instance) :
        """this is the division operator, called upon division of student by number"""
        
        if type(instance) is float or type(instance) is int : # rules for dividing a student by number

            initials = self.initials
            height = self.height / instance 
            birthday = self.birthday / instance
        
        else : # raise error for anything else
            raise Exception('addition not supported for ',type(instance))

        return Student(initials,height,birthday)

    def __floordiv__(self,number) :
        """this is the integer division operator, called upon division of student by int"""
        return self.__div__(number)
    
    def __truediv__(self,number) :
        """this is the float division operator, called upon division of student by float"""
        return self.__div__(number)
    

class StudentArray(ndarray):
    """
    Subclassing of numpy arrays allows us to write additional methods
    on top of the existing library functions. The ndarray.get() method
    allows us to extract an array of attributes from an array of students.
    """

    def __new__(cls, input_array):
        """default initialisation for numpy arrays"""
        obj = np.asarray(input_array).view(cls)
        return obj

    def get(self,attr):
        """custom method to extract any student attributes"""
        return array([ getattr(student,attr) for student in self ])

    def __array_wrap__(self, obj):
        """prevents functions on array from returning zero dimensional objects"""
        if obj.shape == ():
            return obj[()]
        else:
            return np.ndarray.__array_wrap__(self, obj)

    def __array_finalize__(self, obj):
        """numpy needs this. dont worry too much why"""
        if obj is None: return

    
####################### function definitions ############################

def import_students(file_path=None,n_students=5000,mean_height=170.0) :
    """
    returns an array of students given by file, or generates random students.
    
    ---parameters---

    file_path : <str>
        absolute file path to file with three columns separated by a single
        space with initials, heights then birthdays along each one. If none
        provided, then random data is generated.
        
    n_students : <int>
        Number of students to generate if no file is provided.

    """

    if file_path is None : # generate random data

        initials = choice(list(ascii_uppercase),size=n_students)
        heights = normal(loc=mean_height,scale=10.0,size=n_students)
        birthdays = randint(1,366,size=n_students)
    
    else : # import data from file
        initials,heights,birthdays = genfromtxt(file_path,delimiter=' ',dtype=str).T

    return StudentArray([ Student(initial,height,birthday)
                         for initial,height,birthday 
                         in zip(initials,heights,birthdays) ])

def histogram_figure(data,label,**kwargs) :
    
    # compute histogram
    hist,bin_edges = histogram( data, **kwargs )
    bin_centers = (bin_edges[1:]+bin_edges[:-1])/2
    
    # plot it
    plot(bin_centers,hist,label=label)


####################### main program ############################
def main() :
    
    # generate student data, returning StudentArray class
    students = import_students()
    
    # use numpy functions to calculate the mean student
    mean_student = mean(students)
    mean_height = mean_student.height
    mean_birthday = mean_student.birthday
    
    # extract attributes from array subclass
    heights = students.get('height')
    birthdays = students.get('birthday')
    
    # plot results
    figure(figsize=(10,10))
    histogram_figure(heights,'Heights',bins=50,density=True)
    histogram_figure(birthdays,'Birthday',bins=50,density=True)
    
    axvline(mean_height,color='blue')
    axvline(mean_birthday,color='orange')
    
    ylabel(r'Probability Density, $p(x)$',fontsize=16)
    xlabel(r'Bin Center, $x$',fontsize=16)
    
    ylim(0,)
    legend(fontsize=16)
    show()
    
    
# execute main program if this file run directly
if __name__== "__main__":
    main()
