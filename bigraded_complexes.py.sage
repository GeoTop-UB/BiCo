# Bigraded complex
class BigradedComplex():
    def __init__(self, base, dell, delbar, names=None, latex_names=None, CHECK=False):
        self.__base = base
        self.__dimension = {}
        self.__bidegrees = []
        self.__dell = dell
        self.__delbar = delbar
        self.__delldelbar = {}
        self.__names = names
        self.__latex_names = latex_names

        self.__dell_cocycles = {}
        self.__dell_coboundaries = {}
        self.__delbar_cocycles = {}
        self.__delbar_coboundaries = {}
        self.__delldelbar_cocycles = {}
        self.__delldelbar_coboundaries = {}
        self.__dell_and_delbar_cocycles = {}
        self.__dell_and_delbar_coboundaries = {}
        self.__dell_cohomology = {}
        self.__delbar_cohomology = {}
        self.__aeppli_cohomology = {}
        self.__bottchern_cohomology = {}
        self.__reduced_bottchern_cohomology = {}
        self.__reduced_aeppli_cohomology = {}
        self.__zigzags_basis = {}
        self.__squares_basis = {}

        #TODO: Check if the provided bigraded complex is well-defined (i.e. if dell, delbar form a bidifferential) (in case CHECK == True)

        # Record the dimensions of the bigraded components
        for (p,q) in self.__dell:
            if ((p,q) in self.__dimension) == False:
                self.__dimension[(p,q)] = self.__dell[(p,q)].ncols()
        for (p,q) in self.__delbar:
            if ((p,q) in self.__dimension) == False:
                self.__dimension[(p,q)] = self.__delbar[(p,q)].ncols()

        min_p = None
        max_p = None
        min_q = None
        max_q = None
        for (p,q) in self.bidegrees():
            if (min_p == None):
                min_p = p
                max_p = p
                min_q = q
                max_q = q
            else:
                if p < min_p:
                    min_p = p
                if p > max_p:
                    max_p = p
                if q < min_q:
                    min_q = q
                if q > max_q:
                    max_q = q
        self.__min_p = min_p
        self.__max_p = max_p
        self.__min_q = min_q
        self.__max_q = max_q

        # Compute dell-delbar
        for (p,q) in self.__dimension:
            if (p-1, q-1) not in self.__dimension:
                self.__delldelbar[(p-1,q-1)] = Matrix(self.__base, 0, self.__dimension[(p,q)])
            if (p+1,q+1) not in self.__dimension:
                self.__delldelbar[(p,q)] = Matrix(self.__base, 0, self.__dimension[(p,q)])
            elif (p,q+1) not in self.__dimension:
                self.__delldelbar[(p,q)] = Matrix(self.__base, self.__dimension[(p+1,q+1)], self.__dimension[(p,q)])
            else:
                self.__delldelbar[(p,q)] = self.__dell[(p,q+1)] * self.__delbar[(p,q)]

    def base(self):
        r"""
        Return the coefficient field of the `DoubleComplex`.

        EXAMPLES::

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.base()
            Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        return self.__base

    def dimension(self, bidegree=None):
        r"""
        Return the dimension of the ``DoubleComplex`` at the specified bidegree.

        INPUT::

        - ``bidegree`` -- tuple of two integers (default: ``None``)

        OUTPUT::

        Return the dimension of the `DoubleComplex` at bidegree ``bidegree``.
        If ``bidegree`` is set to ``None``, return a dictionary containing
        the dimension of every bigraded component.
        
        EXAMPLES::

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dimension((1,1))
            9
            sage: Iwasawa.dimension()
            {(0, 0): 1,
             (0, 1): 3,
             (0, 2): 3,
             (0, 3): 1,
             (1, 0): 3,
             (1, 1): 9,
             (1, 2): 9,
             (1, 3): 3,
             (2, 0): 3,
             (2, 1): 9,
             (2, 2): 9,
             (2, 3): 3,
             (3, 0): 1,
             (3, 1): 3,
             (3, 2): 3,
             (3, 3): 1}
        """
        if bidegree == None:
            return self.__dimension
        else:
            if bidegree not in self.bidegrees():
                return 0
            else:
                return self.__dimension[bidegree]

    def dell(self, bidegree=None):
        r"""
        Return the matrix of `dell` at the specified bidegree.

        INPUT::
        
        - ``bidegree`` -- tuple of two integers (default: ``None``)

        OUTPUT::

        Return the matrix of `dell` at bidegree ``bidegree``.
        If ``bidegree`` is set to ``None``, return a dictionary containing
        the matrix of `dell` at every bigraded component.
        """
        if bidegree == None:
            return self.__dell
        else:
            if bidegree not in self.bidegrees():
                raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
            else:
                return self.__dell[bidegree]

    def delbar(self, bidegree=None):
        r"""
        Return the matrix of `delbar` at the specified bidegree.

        INPUT::
        
        - ``bidegree`` -- tuple of two integers (default: ``None``)

        OUTPUT::

        Return the matrix of `delbar` at bidegree ``bidegree``.
        If ``bidegree`` is set to ``None``, return a dictionary containing
        the matrix of `delbar` at every bigraded component.
        """
        if bidegree == None:
            return self.__delbar
        else:
            if bidegree not in self.bidegrees():
                raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
            else:
                return self.__delbar[bidegree]

    def delldelbar(self, bidegree=None):
        r"""
        Return the matrix of the composite of `dell delbar` at the specified bidegree.

        INPUT::
        
        - ``bidegree`` -- tuple of two integers (default: ``None``)

        OUTPUT::

        Return the matrix of `dell delbar` at bidegree ``bidegree``.
        If ``bidegree`` is set to ``None``, return a dictionary containing
        the matrix of `dell delbar` at every bigraded component.
        """
        if bidegree == None:
            return self.__delldelbar
        else:
            if bidegree not in self.bidegrees():
                raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
            else:
                return self.__delldelbar[bidegree]

    def names(self, bidegree=None):
        r"""
        Return the names of the canonical basis at the specified bidegree.

        INPUT::
        
        - ``bidegree`` -- tuple of two integers (default: ``None``)

        OUTPUT::

        Return a list with the names of the canonical basis at bidegree ``bidegree``.
        If ``bidegree`` is set to ``None``, return a dictionary containing
        the names of all the canonical bases.

        EXAMPLES:

            sage: KT = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KT.names((1,1))
            [a*abar, b*abar, a*bbar, b*bbar]
            sage: KT.names()
            {(0, 0): [1],
             (0, 1): [abar, bbar],
             (0, 2): [abar*bbar],
             (1, 0): [a, b],
             (1, 1): [a*abar, b*abar, a*bbar, b*bbar],
             (1, 2): [a*abar*bbar, b*abar*bbar],
             (2, 0): [a*b],
             (2, 1): [a*b*abar, a*b*bbar],
             (2, 2): [a*b*abar*bbar]}
        """
        if bidegree == None:
            return self.__names
        else:
            if bidegree not in self.bidegrees():
                raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
            else:
                return self.__names[bidegree]

    def bidegrees(self):
        r"""
        Return a list of all the bidegrees of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.bidegrees()
            [(0, 0),
             (0, 1),
             (0, 2),
             (0, 3),
             (1, 0),
             (1, 1),
             (1, 2),
             (1, 3),
             (2, 0),
             (2, 1),
             (2, 2),
             (2, 3),
             (3, 0),
             (3, 1),
             (3, 2),
             (3, 3)]
        """
        return self.__dimension.keys()

    def zigzags_dimension(self, bidegree=None):
        if bidegree == None:
            dimensions = {}
            for bidegree in self.bidegrees():
                dimensions[bidegree] = self.zigzags_dimension(bidegree=bidegree)
            return dimensions
        else:
            return len(self.zigzags_basis(bidegree))

    def bigraded_component(self, bidegree):
        r"""
        Return the bigraded component with the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.bigraded_component((1,2))
            Free module generated by {a*abar*bbar, b*abar*bbar, c*abar*bbar, a*abar*cbar, b*abar*cbar, c*abar*cbar, a*bbar*cbar, b*bbar*cbar, c*bbar*cbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        if self.names() == None:
            raise AttributeError("The bigraded complex has unspecified names")
        else:
            if bidegree not in self.bidegrees():
                raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
            else:
                return VectorSpace(self.base(), self.names(bidegree))

    def element(self, bidegree, coordinates):
        r"""
        Return the element associated to a bidegree and canonic coordinates.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``coordinates`` -- list of coordinates in the canonical basis of the
        bigraded component of bidegree ``bidegree``

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.element((1,2), (1,0,0,0,1,0,0,-1,0))
            a*abar*bbar + b*abar*cbar - b*bbar*cbar
        """
        if self.names() == None:
            raise AttributeError("The bigraded complex has unspecified names")
        elif bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif len(coordinates) != self.dimension(bidegree):
            raise TypeError("The length of the provided coordinates does not coincide with the dimension of the bigraded component at bidegree " + str(bidegree))
        else:
            return sum(c * b for (c, b) in zip(coordinates, self.names(bidegree)))

################ Dell #################
    def dell_cocycles(self, bidegree, raw=False):
        r"""
        Return the space of `dell` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dell` cocycles of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cocycles((1,2), raw=True)
            Vector space of degree 9 and dimension 6 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
            sage: Iwasawa.dell_cocycles((3,1))
            Free module generated by {a*b*c*abar, a*b*c*bbar, a*b*c*cbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__dell_cocycles:
                self.__dell_cocycles[bidegree] = self.dell(bidegree).right_kernel()
                return self.__dell_cocycles[bidegree]
            else:
                return self.__dell_cocycles[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.dell_cocycles(bidegree, raw=True).basis()])
    
    def dell_cocycles_raw(self, bidegree):
        r"""
        Return the space of raw `dell` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cocycles_raw((1,2))
            Vector space of degree 9 and dimension 6 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        return self.dell_cocycles(bidegree, raw=True)

    def dell_cocycles_inclusion(self, bidegree):
        r"""
        Return an inclusion matrix of the `dell` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cocycles_inclusion((1,3))
            [1 0]
            [0 1]
            [0 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            dell_cocycles = self.dell_cocycles(bidegree, raw=True).basis()
            n = len(dell_cocycles)
            matrix = Matrix(self.base(), self.dimension(bidegree), n)
            for i in range(n):
                for j in range(self.dimension(bidegree)):
                    matrix[j,i] = dell_cocycles[i][j]
            return matrix

    def dell_cocycles_projection(self, bidegree):
        r"""
        Return a projection matrix of the `dell` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cocycles_projection((1,3))
            [1 0 0]
            [0 1 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            dell_cocycles = self.dell_cocycles(bidegree, raw=True)
            basis = list(dell_cocycles.gens())
            n = len(basis)
            for i in range(self.dimension(bidegree)):
                v = vector([int(i == j) for j in range(self.dimension(bidegree))])
                if v not in dell_cocycles:
                    basis += [v]
            change_basis = Matrix(self.base(), self.dimension(bidegree), basis)
            projection = Matrix(self.base(), n, self.dimension(bidegree))
            for i in range(n):
                projection[i,i] = 1
            return projection * change_basis

    def dell_coboundaries(self, bidegree, raw=False):
        r"""
        Return the space of `dell` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dell` coboundaries of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.dell_coboundaries((1,1))
            Free module generated by {a*abar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KodairaThurston.dell_coboundaries((1,1), raw=True)
            Vector space of degree 4 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0]
        """
        (p,q) = bidegree
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__dell_coboundaries:
                if (p-1,q) not in self.bidegrees():
                    self.__dell_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace([0])
                else:
                    self.__dell_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.dell((p-1,q)).transpose(), self.base())
            return self.__dell_coboundaries[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.dell_coboundaries(bidegree, raw=True).basis()])

    def dell_coboundaries_raw(self, bidegree):
        r"""
        Return the space of raw `dell` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_coboundaries_raw((2,2))
            Vector space of degree 9 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]            
        """
        return self.dell_coboundaries(bidegree, raw=True)

############# Delbar ################
    def delbar_cocycles(self, bidegree, raw=False):
        r"""
        Return the space of `delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `delbar` cocycles of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delbar_cocycles((1,3), raw=True)
            Vector space of degree 3 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: Iwasawa.delbar_cocycles((1,1))
            Free module generated by {a*abar, b*abar, c*abar, a*bbar, b*bbar, c*bbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__delbar_cocycles:
                self.__delbar_cocycles[bidegree] = self.delbar(bidegree).right_kernel()
                return self.__delbar_cocycles[bidegree]
            else:
                return self.__delbar_cocycles[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.delbar_cocycles(bidegree, raw=True).basis()])

    def delbar_cocycles_raw(self, bidegree):
        r"""
        Return the space of raw `delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.delbar_cocycles_raw((1,1))
            Vector space of degree 4 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
        """
        return self.delbar_cocycles(bidegree, raw=True)

    def delbar_cocycles_inclusion(self, bidegree):
        r"""
        Return an inclusion matrix of the `delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.delbar_cocycles_inclusion((1,1))
            [1 0 0]
            [0 1 0]
            [0 0 1]
            [0 0 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            delbar_cocycles = self.delbar_cocycles(bidegree, raw=True).basis()
            n = len(delbar_cocycles)
            matrix = Matrix(self.base(), self.dimension(bidegree), n)
            for i in range(n):
                for j in range(self.dimension(bidegree)):
                    matrix[j,i] = delbar_cocycles[i][j]
            return matrix

    def delbar_cocycles_projection(self, bidegree):
        r"""
        Return a projection matrix of the `delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.delbar_cocycles_projection((1,1))
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            delbar_cocycles = self.delbar_cocycles(bidegree, raw=True)
            basis = list(delbar_cocycles.gens())
            n = len(basis)
            for i in range(self.dimension(bidegree)):
                v = vector([int(i == j) for j in range(self.dimension(bidegree))])
                if v not in delbar_cocycles:
                    basis += [v]
            change_basis = Matrix(self.base(), self.dimension(bidegree), basis)
            projection = Matrix(self.base(), n, self.dimension(bidegree))
            for i in range(n):
                projection[i,i] = 1
            return projection * change_basis

    def delbar_coboundaries(self, bidegree, raw=False):
        r"""
        Return the space of `dellbar` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dellbar` coboundaries of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delbar_coboundaries((2,2))
            Free module generated by {a*b*abar*bbar, a*c*abar*bbar, b*c*abar*bbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.delbar_coboundaries((2,2), raw=True)
            Vector space of degree 9 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
        """
        (p,q) = bidegree
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__delbar_coboundaries:
                if (p,q-1) not in self.bidegrees():
                    self.__delbar_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace([0])
                else:
                    self.__delbar_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.delbar((p,q-1)).transpose(), self.base())
            return self.__delbar_coboundaries[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.delbar_coboundaries(bidegree, raw=True).basis()])

    def delbar_coboundaries_raw(self, bidegree):
        r"""
        Return the space of raw `delbar` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delbar_coboundaries_raw((2,2))
            Vector space of degree 9 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
        """
        return self.delbar_coboundaries(bidegree, raw=True)

############### Dell after delbar ################
    def delldelbar_cocycles(self, bidegree, raw=False):
        r"""
        Return the space of `dell delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dell delbar` cocycles of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_cocycles((1,1))
            Free module generated by {a*abar, b*abar, c*abar, a*bbar, b*bbar, c*bbar, a*cbar, b*cbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.delldelbar_cocycles((2,1), raw=True)
            Vector space of degree 9 and dimension 9 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 1]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__delldelbar_cocycles:
                self.__delldelbar_cocycles[bidegree] = self.delldelbar(bidegree).right_kernel()
                return self.__delldelbar_cocycles[bidegree]
            else:
                return self.__delldelbar_cocycles[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.delldelbar_cocycles(bidegree, raw=True).basis()])

    def delldelbar_cocycles_raw(self, bidegree):
        r"""
        Return the space of raw `dell delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_cocycles_raw((1,1))
            Vector space of degree 9 and dimension 8 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        return self.delldelbar_cocycles(bidegree, raw=True)

    def delldelbar_cocycles_inclusion(self, bidegree):
        r"""
        Return an inclusion matrix of the `dell delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_cocycles_inclusion((1,1))
            [1 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0]
            [0 0 0 1 0 0 0 0]
            [0 0 0 0 1 0 0 0]
            [0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            delldelbar_cocycles = self.delldelbar_cocycles(bidegree, raw=True).basis()
            n = len(delldelbar_cocycles)
            matrix = Matrix(self.base(), self.dimension(bidegree), n)
            for i in range(n):
                for j in range(self.dimension(bidegree)):
                    matrix[j,i] = delldelbar_cocycles[i][j]
            return matrix

    def delldelbar_cocycles_projection(self, bidegree):
        r"""
        Return a projection matrix of the `dell delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_cocycles_projection((1,1))
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]        
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            delldelbar_cocycles = self.delldelbar_cocycles(bidegree, raw=True)
            basis = list(delldelbar_cocycles.gens())
            n = len(basis)
            for i in range(self.dimension(bidegree)):
                v = vector([int(i == j) for j in range(self.dimension(bidegree))])
                if v not in delldelbar_cocycles:
                    basis += [v]
            change_basis = Matrix(self.base(), self.dimension(bidegree), basis)
            projection = Matrix(self.base(), n, self.dimension(bidegree))
            for i in range(n):
                projection[i,i] = 1
            return projection * change_basis

    def delldelbar_coboundaries(self, bidegree, raw=False):
        r"""
        Return the space of `dell delbar` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dell delbar` coboundaries of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_coboundaries((2,2))
            Free module generated by {a*b*abar*bbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.delldelbar_coboundaries((2,2), raw=True)
            Vector space of degree 9 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
        """
        (p,q) = bidegree
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__delldelbar_coboundaries:
                if (p-1,q-1) not in self.bidegrees():
                    self.__delldelbar_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace([0])
                else:
                    self.__delldelbar_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.delldelbar((p-1,q-1)).transpose(), self.base())
            return self.__delldelbar_coboundaries[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.delldelbar_coboundaries(bidegree, raw=True).basis()])

    def delldelbar_coboundaries_raw(self, bidegree):
        r"""
        Return the space of raw `dell delbar` coboundaries of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delldelbar_coboundaries_raw((2,2))
            Vector space of degree 9 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
        """
        return self.delldelbar_coboundaries(bidegree, raw=True)

############## Dell intersection delbar ###########
    def dell_and_delbar_cocycles(self, bidegree, raw=False):
        r"""
        Return the intersection of `dell` cocycles with `delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of `dell` cocycles of bidegree ``bidegree`` intersected
        with the vector space of `delbar` cocycles of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.dell_and_delbar_cocycles((1,2))
            Free module generated by {a*abar*bbar, b*abar*bbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KodairaThurston.dell_and_delbar_cocycles((0,2), raw=True)
            Vector space of degree 1 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__dell_and_delbar_cocycles:
                self.__dell_and_delbar_cocycles[bidegree] = self.dell_cocycles(bidegree, raw=True).intersection(self.delbar_cocycles(bidegree, raw=True))
                return self.__dell_and_delbar_cocycles[bidegree]
            else:
                return self.__dell_and_delbar_cocycles[bidegree]
        else:
            return VectorSpace(self.base(), [self.element(bidegree, b) for b in self.dell_and_delbar_cocycles(bidegree, raw=True).basis()])

    def dell_and_delbar_cocycles_raw(self, bidegree):
        r"""
        Return the raw intersection of `dell` and `delbar` cocycles of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_and_delbar_cocycles_raw((1,1))
            Vector space of degree 9 and dimension 4 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]

        """
        return self.dell_and_delbar_cocycles(bidegree, raw=True)

    def dell_and_delbar_cocycles_inclusion(self, bidegree):
        r"""
        Return an inclusion matrix of the intersection of `dell` and `delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_and_delbar_cocycles_inclusion((1,1))
            [1 0 0 0]
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]            
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            dell_and_delbar_cocycles = self.dell_and_delbar_cocycles(bidegree, raw=True).basis()
            n = len(dell_and_delbar_cocycles)
            matrix = Matrix(self.base(), self.dimension(bidegree), n)
            for i in range(n):
                for j in range(self.dimension(bidegree)):
                    matrix[j,i] = dell_and_delbar_cocycles[i][j]
            return matrix

    def dell_and_delbar_cocycles_projection(self, bidegree):
        r"""
        Return a projection matrix of the intersection of `dell` and `delbar` cocycles of the specified bidegree.

        INPUT:

        - bidegree -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_and_delbar_cocycles_projection((1,1))
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]                 
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        else:
            dell_and_delbar_cocycles = self.dell_and_delbar_cocycles(bidegree, raw=True)
            basis = list(dell_and_delbar_cocycles.gens())
            n = len(basis)
            for i in range(self.dimension(bidegree)):
                v = vector([int(i == j) for j in range(self.dimension(bidegree))])
                if v not in dell_and_delbar_cocycles:
                    basis += [v]
            change_basis = Matrix(self.base(), self.dimension(bidegree), basis)
            projection = Matrix(self.base(), n, self.dimension(bidegree))
            for i in range(n):
                projection[i,i] = 1
            return projection * change_basis

############# Cohomologies #############    
    # Cohomology with respect to dell
    def dell_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the cohomology with respect to `dell` at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the cohomology with respect to `dell` of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cohomology_basis((2,3))
            ['[a*c*abar*bbar*cbar]', '[b*c*abar*bbar*cbar]']
            sage: Iwasawa.dell_cohomology_basis((2,3), raw=True)
            [
            (1, 0),
            (0, 1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.dell_cohomology(bidegree, raw=True).basis()
        else:
            dell_cohomology_raw = self.dell_cohomology(bidegree, raw=True)
            raw_basis = [dell_cohomology_raw.lift(b) for b in dell_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def dell_cohomology(self, bidegree, raw=False):
        r"""
        Return the cohomology with respect to `dell` at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Cohomology (as a vector space) with respect to `dell` of bidegree
        ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cohomology((1,3))
            Free module generated by {'[a*abar*bbar*cbar]', '[b*abar*bbar*cbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.dell_cohomology((1,3), raw=True)
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 3 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
            [0 1 0]
            W: Vector space of degree 3 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__dell_cohomology:
                self.__dell_cohomology[bidegree] = self.dell_cocycles(bidegree, raw=True)/self.dell_coboundaries(bidegree, raw=True)
            return self.__dell_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.dell_cohomology_basis(bidegree))

    def dell_cohomology_raw(self, bidegree):
        r"""
        Return the raw cohomology with respect to `dell` of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.dell_cohomology_raw((1,3))
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 3 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
            [0 1 0]
            W: Vector space of degree 3 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        return self.dell_cohomology(bidegree, raw=True)

    # Cohomology with respect to delbar
    def delbar_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the cohomology with respect to `delbar` at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the cohomology with respect to `delbar` of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delbar_cohomology_basis((1,2))
            ['[a*abar*cbar]',
             '[b*abar*cbar]',
             '[c*abar*cbar]',
             '[a*bbar*cbar]',
             '[b*bbar*cbar]',
             '[c*bbar*cbar]']
            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.delbar_cohomology_basis((1,3), raw=True)
            [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.delbar_cohomology(bidegree, raw=True).basis()
        else:
            delbar_cohomology_raw = self.delbar_cohomology(bidegree, raw=True)
            raw_basis = [delbar_cohomology_raw.lift(b) for b in delbar_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def delbar_cohomology(self, bidegree, raw=False):
        r"""
        Return the cohomology with respect to `delbar` at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Cohomology (as a vector space) with respect to `delbar` of bidegree
        ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.delbar_cohomology((1,1))
            Free module generated by {'[b*abar]', '[a*bbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KodairaThurston.delbar_cohomology((1,2), raw=True)
            Vector space quotient V/W of dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 2 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
            [0 1]
            W: Vector space of degree 2 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__delbar_cohomology:
                self.__delbar_cohomology[bidegree] = self.delbar_cocycles(bidegree, raw=True)/self.delbar_coboundaries(bidegree, raw=True)
            return self.__delbar_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.dell_cohomology_basis(bidegree))

    def delbar_cohomology_raw(self, bidegree):
        r"""
        Return the raw cohomology with respect to `delbar` of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.delbar_cohomology_raw((1,1))
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 4 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            W: Vector space of degree 4 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0]
        """
        return self.delbar_cohomology(bidegree, raw=True)

    # Bott-Chern cohomology
    def bottchern_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the Bott-Chern cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the Bott-Chern cohomology of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.bottchern_cohomology_basis((1,1))
            ['[a*abar]', '[b*abar]', '[a*bbar]']
            sage: KodairaThurston.bottchern_cohomology_basis((1,2), raw=True)
            [
            (1, 0),
            (0, 1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.bottchern_cohomology(bidegree, raw=True).basis()
        else:
            bottchern_cohomology_raw = self.bottchern_cohomology(bidegree, raw=True)
            raw_basis = [bottchern_cohomology_raw.lift(b) for b in bottchern_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def bottchern_cohomology(self, bidegree, raw=False):
        r"""
        Return the Bott-Chern cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Bott-Chern cohomology (as a vector space) of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.bottchern_cohomology((1,1))
            Free module generated by {'[a*abar]', '[b*abar]', '[a*bbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KodairaThurston.bottchern_cohomology((1,2), raw=True)
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 2 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
            [0 1]
            W: Vector space of degree 2 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__bottchern_cohomology:
                self.__bottchern_cohomology[bidegree] = (self.dell_cocycles(bidegree, raw=True).intersection(self.delbar_cocycles(bidegree, raw=True)))/self.delldelbar_coboundaries(bidegree, raw=True)
            return self.__bottchern_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.bottchern_cohomology_basis(bidegree))

    def bottchern_cohomology_raw(self, bidegree):
        r"""
        Return the raw Bott-Chern cohomology of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: KodairaThurston = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KodairaThurston.bottchern_cohomology_raw((2,1))
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 2 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
            [0 1]
            W: Vector space of degree 2 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        return self.bottchern_cohomology(bidegree, raw=True)

    # Aeppli cohomology
    def aeppli_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the Aeppli cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the Aeppli cohomology of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.aeppli_cohomology_basis((1,2))
            ['[a*abar*cbar]',
             '[b*abar*cbar]',
             '[c*abar*cbar]',
             '[a*bbar*cbar]',
             '[b*bbar*cbar]',
             '[c*bbar*cbar]']
            sage: Iwasawa.aeppli_cohomology_basis((1,2), raw=True)
            [
            (1, 0, 0, 0, 0, 0),
            (0, 1, 0, 0, 0, 0),
            (0, 0, 1, 0, 0, 0),
            (0, 0, 0, 1, 0, 0),
            (0, 0, 0, 0, 1, 0),
            (0, 0, 0, 0, 0, 1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.aeppli_cohomology(bidegree, raw=True).basis()
        else:
            aeppli_cohomology_raw = self.aeppli_cohomology(bidegree, raw=True)
            raw_basis = [aeppli_cohomology_raw.lift(b) for b in aeppli_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def aeppli_cohomology(self, bidegree, raw=False):
        r"""
        Return the Aeppli cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Aeppli cohomology (as a vector space) of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.aeppli_cohomology((3,1))
            Free module generated by {'[a*b*c*abar]', '[a*b*c*bbar]', '[a*b*c*cbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.aeppli_cohomology((3,2), raw=True)
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 3 and dimension 3 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            W: Vector space of degree 3 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__aeppli_cohomology:
                self.__aeppli_cohomology[bidegree] = self.delldelbar_cocycles(bidegree, raw=True) / (self.dell_coboundaries(bidegree, raw=True) + self.delbar_coboundaries(bidegree, raw=True))
            return self.__aeppli_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.aeppli_cohomology_basis(bidegree))

    def aeppli_cohomology_raw(self, bidegree):
        r"""
        Return the raw Aeppli cohomology of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.aeppli_cohomology_raw((2,2))
            Vector space quotient V/W of dimension 4 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 9 and dimension 9 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 1]
            W: Vector space of degree 9 and dimension 5 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]
        """
        return self.aeppli_cohomology(bidegree, raw=True)

    # Reduced Bott-Chern cohomology
    def reduced_bottchern_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the reduced Bott-Chern cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the reduced Bott-Chern cohomology of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.reduced_bottchern_cohomology_basis((1,2))
            ['[a*abar*bbar]', '[b*abar*bbar]']
            sage: Iwasawa.reduced_bottchern_cohomology_basis((3,2), raw=True)
            [
            (1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.reduced_bottchern_cohomology(bidegree, raw=True).basis()
        else:
            reduced_bottchern_cohomology_raw = self.reduced_bottchern_cohomology(bidegree, raw=True)
            raw_basis = [reduced_bottchern_cohomology_raw.lift(b) for b in reduced_bottchern_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def reduced_bottchern_cohomology(self, bidegree, raw=False):
        r"""
        Return the reduced Bott-Chern cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Reduced Bott-Chern cohomology (as a vector space) of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.reduced_bottchern_cohomology((1,2))
            Free module generated by {'[a*abar*bbar]', '[b*abar*bbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.reduced_bottchern_cohomology((2,3), raw=True)
            Vector space quotient V/W of dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 3 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0]
            W: Vector space of degree 3 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__reduced_bottchern_cohomology:
                self.__reduced_bottchern_cohomology[bidegree] = (self.dell_coboundaries_raw(bidegree).intersection(self.delbar_cocycles_raw(bidegree)) + self.delbar_coboundaries_raw(bidegree).intersection(self.dell_cocycles_raw(bidegree))) / self.delldelbar_coboundaries_raw(bidegree)
            return self.__reduced_bottchern_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.reduced_bottchern_cohomology_basis(bidegree))

    def reduced_bottchern_cohomology_raw(self, bidegree):
        r"""
        Return the raw reduced Bott-Chern cohomology of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: KT = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KT.reduced_bottchern_cohomology_raw((1,1))
            Vector space quotient V/W of dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 4 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0]
            W: Vector space of degree 4 and dimension 0 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            []
        """
        return self.reduced_bottchern_cohomology(bidegree, raw=True)

    # Reduced Aeppli cohomology
    def reduced_aeppli_cohomology_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the reduced Aeppli cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the reduced Aeppli cohomology of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KT = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KT.reduced_aeppli_cohomology_basis((0,1))
            ['[bbar]']
            sage: KT.reduced_aeppli_cohomology_basis((0,1), raw=True)
            [
            (1)
            ]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            return self.reduced_aeppli_cohomology(bidegree, raw=True).basis()
        else:
            reduced_aeppli_cohomology_raw = self.reduced_aeppli_cohomology(bidegree, raw=True)
            raw_basis = [reduced_aeppli_cohomology_raw.lift(b) for b in reduced_aeppli_cohomology_raw.basis()]
            return ['[{}]'.format(self.element(bidegree, b)) for b in raw_basis]

    def reduced_aeppli_cohomology(self, bidegree, raw=False):
        r"""
        Return the reduced Aeppli cohomology at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Reduced Aeppli cohomology (as a vector space) of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:

            sage: KT = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
            sage: KT.reduced_aeppli_cohomology((0,1))
            Free module generated by {'[bbar]'} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KT.reduced_aeppli_cohomology((0,1), raw=True)
            Vector space quotient V/W of dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 2 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
            [0 1]
            W: Vector space of degree 2 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
        """
        if bidegree not in self.bidegrees():
            raise ValueError("The bigraded complex does not have bidegree " + str(bidegree))
        elif raw == True or self.names() == None:
            if bidegree not in self.__reduced_aeppli_cohomology:
                self.__reduced_aeppli_cohomology[bidegree] = self.delldelbar_cocycles(bidegree, raw=True) / (self.dell_coboundaries(bidegree, raw=True) + self.delbar_coboundaries(bidegree, raw=True) + self.dell_cocycles(bidegree, raw=True).intersection(self.delbar_cocycles(bidegree, raw=True)))
            return self.__reduced_aeppli_cohomology[bidegree]
        else:
            return VectorSpace(self.base(), self.reduced_aeppli_cohomology_basis(bidegree))

    def reduced_aeppli_cohomology_raw(self, bidegree):
        r"""
        Return the raw reduced Aeppli cohomology of the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        EXAMPLES:

            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.reduced_aeppli_cohomology_raw((1,2))
            Vector space quotient V/W of dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I where
            V: Vector space of degree 9 and dimension 9 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 1]
            W: Vector space of degree 9 and dimension 7 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        return self.reduced_aeppli_cohomology(bidegree, raw=True)

############ Zig-zags and squares ###############

    # Zigzags
    def zigzags_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the zigzags at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the zigzgas of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:
            
            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.zigzags_basis((2,2))
            [a*c*abar*bbar,
             b*c*abar*bbar,
             a*b*abar*cbar,
             a*c*abar*cbar,
             b*c*abar*cbar,
             a*b*bbar*cbar,
             a*c*bbar*cbar,
             b*c*bbar*cbar]
            sage: Iwasawa.zigzags_basis((1,1), raw=True)
            [(0, 0, 1, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 0, 1, 0, 0, 0),
             (0, 0, 0, 0, 0, 0, 1, 0, 0),
             (0, 0, 0, 0, 0, 0, 0, 1, 0),
             (1, 0, 0, 0, 0, 0, 0, 0, 0),
             (0, 1, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 1, 0, 0, 0, 0, 0),
             (0, 0, 0, 0, 1, 0, 0, 0, 0)]
        """
        if bidegree not in self.__zigzags_basis:
            reduced_aeppli = self.reduced_aeppli_cohomology_raw(bidegree)
            bottchern = self.bottchern_cohomology_raw(bidegree)
            self.__zigzags_basis[bidegree] = [reduced_aeppli.lift(b) for b in reduced_aeppli.basis()] + [bottchern.lift(b) for b in bottchern.basis()]
        if raw == True or self.names() == None:
            return self.__zigzags_basis[bidegree]
        else:
            return [self.element(bidegree, b) for b in self.__zigzags_basis[bidegree]]

    def zigzags(self, bidegree, raw=False):
        r"""
        Return a vector space of zigzags at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of zigzgas of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:
            
            sage: KT.zigzags((1,1))
            Free module generated by {b*bbar, a*abar, b*abar, a*bbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: KT.zigzags((1,2), raw=True)
            Vector space of degree 2 and dimension 2 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [1 0]
            [0 1]
        """
        if raw == True or self.names() == None:
            return VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.zigzags_basis(bidegree, raw=True))
        else:
            return VectorSpace(self.base(), self.zigzags_basis(bidegree))

    # Compute an inclusion of the zigzags into the bigraded component
    def zigzags_inclusion(self, bidegree):
        zigzags = self.zigzags_basis(bidegree, raw=True)
        return Matrix(self.__base, len(zigzags), self.__dimension[bidegree], [z for z in zigzags]).transpose()

    # Compute the projection form the bigraded component to the zigzags
    def zigzags_projection(self, bidegree):
        zigzags = self.zigzags_basis(bidegree, raw=True)
        n = len(zigzags)
        change_basis = Matrix(self.base(), self.dimension(bidegree), [z for z in zigzags] + [s for s in self.squares_basis(bidegree, raw=True)]).transpose().inverse()
        projection = Matrix(self.base(), n, self.dimension(bidegree))
        for i in range(n):
            projection[i,i] = 1
        return projection*change_basis

    # Return the deformation retract associated to the zigzags
    # Output: zigzags as a bigraded complex (BigradedComplex), inclusion (BigradedComplexMap), projection (BigradedComplexMap), homotopy (BigradedComplexMap)
    def zigzags_deformation_retract(self):
        names = {}
        dell = {}
        delbar = {}
        inclusion = {}
        projection = {}
        homotopy = {}

        for (p,q) in self.bidegrees():
            names[(p,q)] = self.zigzags_basis((p,q))
            if (p+1,q) in self.bidegrees():
                dell[(p,q)] = self.zigzags_projection((p+1,q)) * self.dell((p,q)) * self.zigzags_inclusion((p,q))
            else:
                dell[(p,q)] = Matrix(self.base(), 0, len(names[(p,q)]))
            if (p,q+1) in self.bidegrees():
                delbar[(p,q)] = self.zigzags_projection((p,q+1)) * self.delbar((p,q)) * self.zigzags_inclusion((p,q))
            else:
                delbar[(p,q)] = Matrix(self.base(), 0, len(names[(p,q)]))
            inclusion[(p,q)] = self.zigzags_inclusion((p,q))
            projection[(p,q)] = self.zigzags_projection((p,q))
            if (p-1,q-1) in self.bidegrees():
                image_delldelbar = self.delldelbar((p-1,q-1)).columns()
                basis = []
                preimage = []
                for i in range(self.dimension((p-1,q-1))):
                    if image_delldelbar[i] not in VectorSpace(self.base(), self.dimension((p,q))).subspace(basis):
                        basis += [image_delldelbar[i]]
                        preimage += [vector([int(j == i) for j in range(self.dimension((p-1,q-1)))])]
                for b in self.squares_basis((p,q), raw=True):
                    if b not in VectorSpace(self.base(), self.dimension((p,q))).subspace(basis):
                        basis += [b]
                for i in range(self.dimension((p,q))):
                    v = vector((int(i == j) for j in range(self.dimension((p,q)))))
                    if v not in VectorSpace(self.base(), self.dimension((p,q))).subspace(basis):
                        basis += [v]
                change_basis = Matrix(self.base(), self.dimension((p,q)), basis).transpose().inverse()
                h = Matrix(self.base(), self.dimension((p,q)), self.dimension((p-1,q-1)))
                for i in range(len(preimage)):
                    h[i,:] = preimage[i]
                h = h.transpose()
                homotopy[(p,q)] = h*change_basis

            else:
                homotopy[(p,q)] = Matrix(self.base(), self.dimension((p,q)), 0)

        zigzags_bicpx = BigradedComplex(self.base(), dell, delbar, names=names)

        return zigzags_bicpx, BigradedComplexMap(zigzags_bicpx, self, inclusion), BigradedComplexMap(self, zigzags_bicpx, projection), BigradedComplexMap(self, self, homotopy, bidegree=(-1,-1))

    # Squares
    def squares_basis(self, bidegree, raw=False):
        r"""
        Return a basis for the squares at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Basis for the squares of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis is made of coordinate vectors. Otherwise, the basis
        depends on the names of the bigraded complex.

        EXAMPLES:
            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.squares_basis((2,1))
            [a*b*cbar]
            sage: Iwasawa.squares_basis((2,2), raw=True)
            [(1, 0, 0, 0, 0, 0, 0, 0, 0)]
        """
        if bidegree not in self.__squares_basis:
            zigzags = self.zigzags(bidegree, raw=True)
            self.__squares_basis[bidegree] = []
            for i in range(self.dimension(bidegree)):
                v = vector([int(i==j) for j in range(self.dimension(bidegree))])
                if v not in zigzags:
                    self.__squares_basis[bidegree] += [v]
        if raw == True or self.names() == None:
            return self.__squares_basis[bidegree]
        else:
            return [self.element(bidegree, b) for b in self.squares_basis(bidegree, raw=True)]

    def squares(self, bidegree, raw=False):
        r"""
        Return a vector space of squares at the specified bidegree.

        INPUT:

        - ``bidegree`` -- tuple of two integers

        - ``raw`` -- boolean (default: ``False``)

        OUTPUT:

        - Vector space of squares of bidegree ``bidegree``.
        If ``raw`` is set to ``True`` (or the bigraded complex has unspecified
        names), the basis of the vector space is made of coordinate vectors.
        Otherwise, the basis depends on the names of the bigraded complex.

        EXAMPLES:
            sage: Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
            sage: Iwasawa.squares((2,1))
            Free module generated by {a*b*cbar} over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            sage: Iwasawa.squares((1,2), raw=True)
            Vector space of degree 9 and dimension 1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
            Basis matrix:
            [0 0 1 0 0 0 0 0 0]
        """
        if raw == True or self.__names == None:
            return VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.squares_basis(bidegree, raw=True))
        else:
            return VectorSpace(self.base(), self.squares_basis(bidegree))

############# Ascii art ################

    # Method used to produce the asii art tables for the cohomologies
    # Data is a dictionary where each (p,q)-entry contains a basis for the (p,q)-th component
    # n_generators_row is an integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1, there will be no limit
    def __ascii_art_table(self, data, n_generators_row):
        width = {}
        height = {}
        total_width = 0
        string = {}

        for (p,q) in self.__dimension:
            n = len(data[(p,q)])
            string[(p,q)] = [""]
            jump = 0
            for i in range(n-1):
                string[(p,q)][jump] = string[(p,q)][jump] + str(data[(p,q)][i]) + ", "
                if i%n_generators_row == n_generators_row-1:
                    string[(p,q)][jump] = string[(p,q)][jump][:-1]
                    jump += 1
                    string[(p,q)] += [""]
            if n != 0:
                string[(p,q)][jump] = string[(p,q)][jump] + str(data[(p,q)][n-1])
            if (q in height) == False:
                height[q] = jump+1
            elif height[q] < jump+1:
                height[q] = jump+1

        for (p,q) in self.bidegrees():
            for h in range(len(string[(p,q)])):
                if (p in width) == False:
                    width[p] = len(string[(p,q)][h])
                elif width[p] < len(string[(p,q)][h]):
                    width[p] = len(string[(p,q)][h])

        number_length = len(str(self.__min_p))
        if len(str(self.__max_p)) > number_length:
            number_length = len(str(self.__max_p))
        if len(str(self.__min_q)) > number_length:
            number_length = len(str(self.__min_q))
        if len(str(self.__max_q)) > number_length:
            number_length = len(str(self.__max_q))

        for p in range(self.__min_p, self.__max_p+1):
            if (p in width) == False:
                width[p] = number_length
            elif width[p] < number_length:
                width[p] = number_length
            total_width += width[p]

        for q in range(self.__min_q, self.__max_q+1):
            if (q in height) == False:
                height[q] = 1

        # Horizontal line
        hline = ["-"]*(number_length + total_width + 5*(self.__max_p - self.__min_p + 1) + 2)
        index = 2+number_length
        for p in range(self.__min_p, self.__max_p+1):
            hline[index] = "+"
            index += width[p] + 5
        hline = "".join(hline)

        for minus_q in range(self.__max_q-self.__min_q+1):
            q = self.__max_q-minus_q
            line = [""]*height[q]
            for h in range(height[q]):
                if (h == 0):
                    line[h] = str(q) + " "*(number_length - len(str(q)))
                else:
                    line[h] = " "*number_length
                for p in range(self.__min_p,  self.__max_p+1):
                    if (p,q) in string:
                        if h < len(string[(p,q)]):
                            line[h] += "  |  " + string[(p,q)][h] + " "*(width[p] - len(string[(p,q)][h]))
                        else:
                            line[h] += "  |  " + " "*width[p]
                    else:
                        line[h] += "  |  " + " "*width[p]
                print(line[h])
            print(hline)

        line = " "*number_length
        for p in range(self.__min_p, self.__max_p+1):
            line += "  |  " + str(p) + " "*(width[p] - len(str(p)))
        print(line)

    # Dell cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_dell_cohomology(self, n_generators_row = -1):
        print("Anti-Dolbeault cohomology:\n")
        self.__ascii_art_table({bidegree: self.dell_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)
        

    # Delbar cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_delbar_cohomology(self, n_generators_row = -1):
        print("Dolbeault cohomology:\n")
        self.__ascii_art_table({bidegree: self.delbar_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

    # Bott-Chern cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_bottchern_cohomology(self, n_generators_row = -1):
        print("Bott-Chern cohomology:\n")
        self.__ascii_art_table({bidegree: self.bottchern_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)
        

    # Aeppli cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_aeppli_cohomology(self, n_generators_row = -1):
        print("Aeppli cohomology:\n")
        self.__ascii_art_table({bidegree: self.aeppli_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

    # Reduced Bott-Chern cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_reduced_bottchern_cohomology(self, n_generators_row = -1):
        print("Reduced Bott-Chern cohomology:\n")
        self.__ascii_art_table({bidegree: self.reduced_bottchern_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

    # Reduced Aeppli cohomology ascii art
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_reduced_aeppli_cohomology(self, n_generators_row = -1):
        print("Reduced Aeppli cohomology:\n")
        self.__ascii_art_table({bidegree: self.reduced_aeppli_cohomology_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

    # Displays ascii art for Dolbeault, anti-Dolbeault, Bott-Chern and Aeppli cohomologies
    # n_generators_row is an optional integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1 (default), there will be no limit
    def _ascii_art_cohomologies(self, n_generators_row = -1):
        self._ascii_art_delbar_cohomology(n_generators_row)
        print("\n")
        self._ascii_art_dell_cohomology(n_generators_row)
        print("\n")
        self._ascii_art_bottchern_cohomology(n_generators_row)
        print("\n")
        self._ascii_art_aeppli_cohomology(n_generators_row)

    def _ascii_art_zigzags(self, n_generators_row = -1):
        print("Zig-zags:\n")
        self.__ascii_art_table({bidegree: self.zigzags_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

    def _ascii_art_squares(self, n_generators_row = -1):
        print("Squares:\n")
        self.__ascii_art_table({bidegree: self.squares_basis(bidegree) for bidegree in self.bidegrees()}, n_generators_row)

# BigradedComplexMap
class BigradedComplexMap():
    def __init__(self, domain, codomain, mapp, bidegree=(0,0)):
        self.__domain = domain
        self.__codomain = codomain
        self.__base = domain.base()
        self.__map = mapp
        self.__bidegree = bidegree
        self.__dell_boundary = {}
        self.__delbar_boundary = {}
        self.__induced_dell_cohomology = {}
        self.__induced_delbar_cohomology = {}
        self.__induced_bottchern_cohomology = {}
        self.__induced_aeppli_cohomology = {}
        self.__is_morphism = None
        self.__kernel = None

    # Call
    def __call__(self, bidegree, element, raw=False):
        if bidegree not in self.__domain.bidegrees():
            raise ValueError("The input bidegree " + str(bidegree) + " must be a bidegree in the domain.")
        (r,s) = self.__bidegree
        (p,q) = bidegree
        if (p+r,q+s) in self.__codomain.bidegrees():
            if raw == True or self.__codomain.names() == None:
                return self.__map[bidegree] * vector(element)
            else:
                return self.__codomain.element((p+r,q+s), self.__call__(bidegree, element, raw=True))
        else:
            return 0

    #Various pieces of information
    def domain(self):
        return self.__domain
    def codomain(self):
        return self.__codomain
    def bidegree(self):
        return self.__bidegree
    def base(self):
        return self.__base
    def map(self, bidegree):
        return self.__map[bidegree]

    # Compute the dell-boundary of the map, i.e. [dell, f]
    def dell_boundary(self, bidegree = None):
        if bidegree == None:
            for bidegree in self.__map:
                if (bidegree in self.__dell_boundary) == False:
                    self.__dell_boundary[bidegree] = self.dell_boundary(bidegree)
            return self.__dell_boundary
        else:
            (p,q) = bidegree
            (r,s) = self.__bidegree
            if (bidegree in self.__dell_boundary) == False:
                if (p+r+1,q+s) in self.__codomain.bidegrees():
                    self.__dell_boundary[bidegree] = Matrix(self.__base, self.__codomain.dimension((p+r+1,q+s)), self.__domain.dimension(bidegree))
                    if (p+r,q+s) in self.__codomain.bidegrees():
                        self.__dell_boundary[bidegree] += self.__codomain.dell((p+r,q+s)) * self.__map[bidegree]
                    if (p+1,q) in self.__domain.bidegrees():
                        self.__dell_boundary[bidegree] += (-1+2*((self.__bidegree[0] + self.__bidegree[1])%2))*self.__map[(p+1,q)] * self.__domain.dell(bidegree)
                else:
                    self.__dell_boundary[bidegree] = Matrix(self.__base, 0, self.__domain.dimension(bidegree))
            return self.__dell_boundary[bidegree]

    # Compute the delbar-boundary of the map, i.e. [delbar, f]
    def delbar_boundary(self, bidegree = None):
        if bidegree == None:
            for bidegree in self.__map:
                if (bidegree in self.__delbar_boundary) == False:
                    self.__delbar_boundary[bidegree] = self.delbar_boundary(bidegree)
            return self.__delbar_boundary
        else:
            (p,q) = bidegree
            (r,s) = self.__bidegree
            if (bidegree in self.__delbar_boundary) == False:
                if (p+r,q+s+1) in self.__codomain.bidegrees():
                    self.__delbar_boundary[bidegree] = Matrix(self.__base, self.__codomain.dimension((p+r,q+s+1)), self.__domain.dimension(bidegree))
                    if (p+r,q+s) in self.__codomain.bidegrees():
                        self.__delbar_boundary[bidegree] += self.__codomain.delbar((p+r,q+s)) * self.__map[bidegree]
                    if (p,q+1) in self.__domain.bidegrees():
                        self.__delbar_boundary[bidegree] += (-1+2*((self.__bidegree[0] + self.__bidegree[1])%2))*self.__map[(p,q+1)] * self.__domain.delbar(bidegree)
                else:
                    self.__delbar_boundary[bidegree] = Matrix(self.__base, 0, self.__domain.dimension(bidegree))
            return self.__delbar_boundary[bidegree]

    # Check if the map is a morphism
    def is_morphism(self):
        if self.__is_morphism == None:
            if self.__bidegree == (0,0):
                self.__is_morphism = True
                dell_boundary = self.dell_boundary()
                delbar_boundary = self.delbar_boundary()
                for (p,q) in self.__map:
                    if ((p+1, q) in self.__codomain.dimension()) == False:
                        if dell_boundary[(p,q)] != Matrix(self.__base, 0, self.__domain.dimension((p,q))):
                            self.__is_morphism = False
                    elif dell_boundary[(p,q)] != Matrix(self.__base, self.__codomain.dimension((p+1,q)), self.__domain.dimension((p,q))):
                        self.__is_morphism = False
                    if ((p, q+1) in self.__codomain.dimension()) == False:
                        if delbar_boundary[(p,q)] != Matrix(self.__base, 0, self.__domain.dimension((p,q))):
                            self.__is_morphism = False
                    elif delbar_boundary[(p,q)] != Matrix(self.__base, self.__codomain.dimension((p,q+1)), self.__domain.dimension((p,q))):
                        self.__is_morphism = False
            else:
                self.__is_morphism = False
        return self.__is_morphism

    # ascii art representation of the map
    def _ascii_art_map(self):
        max_length = 0
        for bidegree in self.__domain.names():
            for name in self.__domain.names()[bidegree]:
                if len(str(name)) > max_length:
                    max_length = len(str(name))

        (r,s) = self.__bidegree
        for bidegree in self.__domain.dimension():
            (p,q) = bidegree
            has_columns = (p+r,q+s) in self.__codomain.bidegrees()
            print("Bidegree " + str(bidegree) + ":")
            for i in range(self.__domain.dimension(bidegree)):
                x = self.__domain.names()[bidegree][i]
                length = len(str(x))
                if has_columns:
                    print("\t" + str(x) + " "*(max_length+2-length) + "|---->  " + str(self.__codomain.element((p+r,q+s), self.__map[bidegree].column(i))))
                else:
                    print("\t" + str(x) + " "*(max_length+2-length) + "|---->  0")
            print("\n")

# BidifferentialBigradedAlgebra
class BidifferentialBigradedCommutativeAlgebra(BigradedComplex):
    def __init__(self, algebra, dell_dictionary, delbar_dictionary, min_deg, max_deg):    
        self.__zigzags = None
        self.__zigzags_i = None
        self.__zigzags_p = None
        self.__zigzags_h11 = None
        self.__zigzags_h10 = None
        self.__zigzags_h01 = None
        self.__algebra = algebra
        self.__min_deg = min_deg
        self.__max_deg = max_deg

        dictionary_is_matrix = (dell_dictionary[list(dell_dictionary.keys())[0]] not in algebra)

        bigraded_basis = {}
        if dictionary_is_matrix:
            dell_matrix = dell_dictionary
            delbar_matrix = delbar_dictionary
        else:
            dell_matrix = {}
            delbar_matrix = {}
        for p in range(self.__min_deg, self.__max_deg+1):
            for q in range(self.__min_deg, self.__max_deg+1):
                basis = self.__algebra.basis((p,q))
                if basis != []:
                    bigraded_basis[(p,q)] = basis

                    # Compute the differential matrices
                    if not dictionary_is_matrix:
                        dell_matrix[(p,q)] = self.__algebra.differential(dell_dictionary).differential_matrix_multigraded((p,q)).transpose()
                        delbar_matrix[(p,q)] = self.__algebra.differential(delbar_dictionary).differential_matrix_multigraded((p,q)).transpose()

        BigradedComplex.__init__(self, self.__algebra.base(), dell_matrix, delbar_matrix, names=bigraded_basis)

    def algebra(self):
        return self.__algebra
    def min_deg(self):
        return self.__min_deg
    def max_deg(self):
        return self.__max_deg
    def zigzags_i(self):
        if self.__zigzags_i == None:
            self.__compute_zigzags_deformation_retract()
        return self.__zigzags_i
    def zigzags_p(self):
        if self.__zigzags_p == None:
            self.__compute_zigzags_deformation_retract()
        return self.__zigzags_p
    def zigzags_h11(self):
        if self.__zigzags_h11 == None:
            self.__compute_zigzags_deformation_retract()
        return self.__zigzags_h11
    def zigzags_h10(self):
        if self.__zigzags_h10 == None:
            self.__compute_zigzags_deformation_retract()
        return self.__zigzags_h10
    def zigzags_h01(self):
        if self.__zigzags_h01 == None:
            self.__compute_zigzags_deformation_retract()
        return self.__zigzags_h01

    def subalgebra(self, subalgebra_generators):
        # Compute the differentials of the generators (which will also become generators), and their bidegrees
        # Write the names of the generators (they will be used as their names in the subalgebra)
        dell_generators = []
        delbar_generators = []
        delldelbar_generators = []                
        names = ["x" + str(i) for i in range(len(subalgebra_generators))]
        dell_names = []
        delbar_names = []
        delldelbar_names = []
        bidegrees = [generator.degree() for generator in subalgebra_generators]
        dell_bidegrees = []
        delbar_bidegrees = []
        delldelbar_bidegrees = []
        index = 0
        for generator in subalgebra_generators:
            (p,q) = generator.degree()
            if (p+1,q) in self.bidegrees() and self.dell((p,q))*vector(generator.basis_coefficients()) != 0:
                dell_generators += [self.element((p+1,q), self.dell((p,q))*vector(generator.basis_coefficients()))]
                dell_names += ["y" + str(index)]
                dell_bidegrees += [(p+1,q)]
            if (p,q+1) in self.bidegrees() and self.delbar((p,q))*vector(generator.basis_coefficients()) != 0:
                delbar_generators += [self.element((p,q+1), self.delbar((p,q))*vector(generator.basis_coefficients()))]
                delbar_names += ["z" + str(index)]
                delbar_bidegrees += [(p,q+1)]
            if (p+1,q+1) in self.bidegrees() and self.delldelbar((p,q))*vector(generator.basis_coefficients()) != 0:
                delldelbar_generators += [self.element((p+1,q+1), self.delldelbar((p,q))*vector(generator.basis_coefficients()))]
                delldelbar_names += ["w" + str(index)]
                delldelbar_bidegrees += [(p+1,q+1)]
            index += 1
        bidegrees = bidegrees + dell_bidegrees + delbar_bidegrees + delldelbar_bidegrees
        generators = subalgebra_generators + dell_generators + delbar_generators + delldelbar_generators
        names = names + dell_names + delbar_names + delldelbar_names

        # Order the generators according to their total degree (required for GradedCommutativeAlgebra)
        total_degrees = {}
        for i in range(len(generators)):
            (p,q) = bidegrees[i]
            if p+q in total_degrees:
                total_degrees[p+q] += [i]
            else:
                total_degrees[p+q] = [i]
        ordered_total_degrees = list(total_degrees.keys())
        ordered_total_degrees.sort()
        ordered_bidegrees = []
        ordered_generators = []
        names_str = ""
        for total_degree in ordered_total_degrees:
            for i in total_degrees[total_degree]:
                ordered_bidegrees += [generators[i].degree()]
                ordered_generators += [generators[i]]
                names_str += names[i] + ","
        names_str = names_str[:-1]

        # Algebra that maps to the initial algebra. The quotient of A by the kernel of the map will give the subalgebra
        A = GradedCommutativeAlgebra(self.__algebra.base(), names=names_str, degrees=ordered_bidegrees)

        # Compute the kernel of the map A --> algebra
        inclusion = Hom(A, self.algebra())(ordered_generators)
        inclusion_matrix = {}
        section_matrix = {}
        ideal = []
        for bidegree in self.bidegrees():
            image_space = VectorSpace(self.base(), self.algebra().basis(bidegree))
            inclusion_matrix[bidegree] = Matrix(self.base(), [inclusion(u).basis_coefficients() if inclusion(u) != 0 else vector([0 for _ in range(self.dimension(bidegree))]) for u in A.basis(bidegree)]).transpose()
            section_matrix[bidegree] = inclusion_matrix[bidegree].pseudoinverse()
            kernel = inclusion_matrix[bidegree].right_kernel().basis()
            ideal += [sum(c*b for (c,b) in zip(u, A.basis(bidegree))) for u in kernel]

        # Quotient
        B = A.quotient(A.ideal(ideal))
        projection = {}
        lift = {}
        for bidegree in self.bidegrees():
            if A.basis(bidegree) != []:
                if B.basis(bidegree) == []:
                    projection[bidegree] = Matrix(self.base(), 0, len(A.basis(bidegree)))
                else:
                    projection[bidegree] = Matrix(self.base(), [B(u).basis_coefficients() if B(u) != 0 else vector(0 for _ in range(len(B.basis(bidegree)))) for u in A.basis(bidegree)]).transpose()
            if B.basis(bidegree) != []:
                if A.basis(bidegree) == []:
                    lift[bidegree] = Matrix(self.base(), 0, len(B.basis(bidegree)))
                else:
                    lift[bidegree] = Matrix(self.base(), [A(u.lift()).basis_coefficients() if u.lift() != 0 else vector([0 for _ in range(len(A.basis(bidegree)))]) for u in B.basis(bidegree)]).transpose()

        # Compute dell and delbar
        dell = {}
        delbar = {}
        for (p,q) in self.bidegrees():
            if B.basis((p,q)) != []:
                if B.basis((p+1,q)) == []:
                    dell[(p,q)] = Matrix(self.base(), 0, len(B.basis((p,q))))
                else:
                    dell[(p,q)] = projection[(p+1,q)]*section_matrix[(p+1,q)]*self.dell((p,q))*inclusion_matrix[(p,q)]*lift[(p,q)]
                if B.basis((p,q+1)) == []:
                    delbar[(p,q)] = Matrix(self.base(), 0, len(B.basis((p,q))))
                else:
                    delbar[(p,q)] = projection[(p,q+1)]*section_matrix[(p,q+1)]*self.delbar((p,q))*inclusion_matrix[(p,q)]*lift[(p,q)]

        # Dictionary that realizes the elements in B as elements in the original algebra
        dictionary = {}
        for bidegree in self.bidegrees():
            for u in B.basis(bidegree):
                dictionary[u] = inclusion(A(u.lift()))

        return BidifferentialBigradedCommutativeAlgebra(B, dell, delbar, self.__min_deg, self.__max_deg), dictionary

    # Produces a BidifferentialBigradedCommutativeAlgebra from a nilmanifold and an almost complex structure
    # Caution! The almost complex structure is assumed to be integrable
    @staticmethod
    def from_nilmanifold(lie_algebra, ac_structure, labels=None, normalization_coefficients=None, latex_generators=None):
        original_generators = lie_algebra.gens()
        dimension = len(original_generators)
        if normalization_coefficients == None:
            normalization_coefficients = [1 for _ in range(dimension)]
        if dimension % 2 == 1:
            raise "The Lie algebra must be even-dimensional."
        elif dimension == 0:
            return BidifferentialBigradedCommutativeAlgebra.unit(lie_algebra.base())
        else:
            if labels == None:
                labels = ['a%s' %j for j in range(dimension/2)]+['b%s' %j for j in range(dimension/2)]

            eigenvectors = ac_structure.eigenvectors_right()
            if (eigenvectors[0][0] == I):
                basis_coefficients = eigenvectors[0][1] + eigenvectors[1][1]
            else:
                basis_coefficients = eigenvectors[1][1] + eigenvectors[0][1]
            
            basis = [normalization_coefficients[i]*sum(basis_coefficients[i][j] * original_generators[j] for j in range(dimension)) for i in range(dimension)]
            change_basis = Matrix(lie_algebra.base(), dimension, dimension, basis_coefficients).transpose().inverse()

            # Define GradedCommutativeAlgebra
            algebra = GradedCommutativeAlgebra(lie_algebra.base(), names=labels, degrees=tuple([(1,0) for _ in range(dimension/2)]+[(0,1) for _ in range(dimension/2)]))
            generators = algebra.gens()

            # Compute the differentials of the Chevalley-Eilenberg bigraded algebra
            dell_dict = {g: 0 for g in generators}
            delbar_dict = {g: 0 for g in generators}
            value = {}
            for k in range(dimension/2):
                for i in range(dimension/2):
                    for j in range(i+1, dimension/2):
                        if (i,j) not in value:
                            value[(i,j)] = -change_basis*vector(lie_algebra.bracket(basis[i], basis[j]))
                        dell_dict[generators[k]] += value[(i,j)][k]*generators[i]*generators[j]
                    for j in range(dimension/2, dimension):
                        if (i,j) not in value:
                            value[(i,j)] = -change_basis*vector(lie_algebra.bracket(basis[i], basis[j]))
                        delbar_dict[generators[k]] += value[(i,j)][k]*generators[i]*generators[j]
            for k in range(dimension/2, dimension):
                for i in range(dimension/2):
                    for j in range(dimension/2, dimension):
                        if (i,j) not in value:
                            value[(i,j)] = -change_basis*vector(lie_algebra.bracket(basis[i], basis[j]))
                        dell_dict[generators[k]] += value[(i,j)][k]*generators[i]*generators[j]
                for i in range(dimension/2, dimension):
                    for j in range(i+1, dimension):
                        if (i,j) not in value:
                            value[(i,j)] = -change_basis*vector(lie_algebra.bracket(basis[i], basis[j]))
                        delbar_dict[generators[k]] += value[(i,j)][k]*generators[i]*generators[j]

            return BidifferentialBigradedCommutativeAlgebra(algebra, dell_dict, delbar_dict, 0, dimension)

    def __compute_zigzags_deformation_retract(self):
        self.__zigzags, self.__zigzags_i, self.__zigzags_p, self.__zigzags_h11 = self.zigzags_deformation_retract()
        self.__zigzags_h10 = BigradedComplexMap(self, self, self.__zigzags_h11.delbar_boundary(), bidegree=(-1,0))
        self.__zigzags_h01 = self.__zigzags_h11.dell_boundary()
        for bidegree in self.__zigzags_h01:
            self.__zigzags_h01[bidegree] = -self.__zigzags_h01[bidegree]
        self.__zigzags_h01 = BigradedComplexMap(self, self, self.__zigzags_h01, bidegree=(0,-1))

    # Returns the dictionary of the specified operation. If a tuple of elements does not appear as a key in the dictionary, the operation vanishes on it
    # If compute_symmetries is set to False, the superfluous operations will not be computed, knowing that: m3^{-1,-1}(x,y,z) = +/-\overline{m3^{-1,-1}(\bar x, \bar y, \bar z)} and m3^{-1,-1}(x,y,z) = -m3^{-1,-1}(z,y,x).
    def operation(self, arity, operation_bidegree, show_progress=False, compute_symmetries=True):
        if self.__zigzags == None:
            self.__compute_zigzags_deformation_retract()
        dictionary = {}
        total_bidegrees = Tuples(self.bidegrees(), arity)

        if show_progress == True:
            n = 0
            n_total_bidegrees = len(total_bidegrees)
            previous_progress = 0

        if arity == 3 and operation_bidegree == (-1,0):
            lifted_product = {}
            for bidegrees in total_bidegrees:
                if compute_symmetries or bidegrees[0][0] + bidegrees[1][0] + bidegrees[2][0] <= bidegrees[0][1] + bidegrees[1][1] + bidegrees[2][1]:
                    for elements in cartesian_product([self.zigzags_basis(bidegree) for bidegree in bidegrees]):
                        result_bidegree = (bidegrees[0][0]+bidegrees[1][0]+bidegrees[2][0]-1, bidegrees[0][1]+bidegrees[1][1]+bidegrees[2][1])
                        if result_bidegree in self.bidegrees() and (compute_symmetries or (elements[2], elements[1], elements[0]) not in dictionary):
                            if (elements[0], elements[1]) not in lifted_product:
                                prod = elements[0]*elements[1]
                                if prod == 0 or (bidegrees[0][0]+bidegrees[1][0], bidegrees[0][1]+bidegrees[1][1]) not in self.bidegrees():
                                    lifted_product[(elements[0], elements[1])] = 0
                                else:
                                    lifted_product[(elements[0], elements[1])] = self.__zigzags_h10((bidegrees[0][0] + bidegrees[1][0], bidegrees[0][1] + bidegrees[1][1]), (elements[0]*elements[1]).basis_coefficients())
                            if (elements[1], elements[2]) not in lifted_product:
                                prod = elements[1]*elements[2]
                                if prod == 0 or (bidegrees[1][0]+bidegrees[2][0], bidegrees[1][1]+bidegrees[2][1]) not in self.bidegrees():
                                    lifted_product[(elements[1], elements[2])] = 0
                                else:
                                    lifted_product[(elements[1], elements[2])] = self.__zigzags_h10((bidegrees[1][0] + bidegrees[2][0], bidegrees[1][1] + bidegrees[2][1]), (elements[1]*elements[2]).basis_coefficients())
                            result = lifted_product[(elements[0], elements[1])]*elements[2] - elements[0]*lifted_product[(elements[1], elements[2])]
                            if result != 0:
                                projected_result = self.__zigzags_p(result_bidegree, result.basis_coefficients())
                                if projected_result != 0:                            
                                    dictionary[tuple(elements)] = projected_result

                if show_progress == True:
                    n += 1
                    progress = int(100*n/n_total_bidegrees)
                    if progress > previous_progress:
                        print("Progress: " + str(progress) + "%")
                    previous_progress = progress

        elif arity == 3 and operation_bidegree == (0,-1):
            lifted_product = {}
            for bidegrees in total_bidegrees:
                if compute_symmetries or bidegrees[0][0] + bidegrees[1][0] + bidegrees[2][0] <= bidegrees[0][1] + bidegrees[1][1] + bidegrees[2][1]:
                    for elements in cartesian_product([self.zigzags_basis(bidegree) for bidegree in bidegrees]):
                        result_bidegree = (bidegrees[0][0]+bidegrees[1][0]+bidegrees[2][0], bidegrees[0][1]+bidegrees[1][1]+bidegrees[2][1]-1)
                        if result_bidegree in self.bidegrees() and (compute_symmetries or (elements[2], elements[1], elements[0]) not in dictionary):
                            if (elements[0], elements[1]) not in lifted_product:
                                prod = elements[0]*elements[1]
                                if prod == 0 or (bidegrees[0][0]+bidegrees[1][0], bidegrees[0][1]+bidegrees[1][1]) not in self.bidegrees():
                                    lifted_product[(elements[0], elements[1])] = 0
                                else:
                                    lifted_product[(elements[0], elements[1])] = self.__zigzags_h01((bidegrees[0][0] + bidegrees[1][0], bidegrees[0][1] + bidegrees[1][1]), (elements[0]*elements[1]).basis_coefficients())
                            if (elements[1], elements[2]) not in lifted_product:
                                prod = elements[1]*elements[2]
                                if prod == 0 or (bidegrees[1][0]+bidegrees[2][0], bidegrees[1][1]+bidegrees[2][1]) not in self.bidegrees():
                                    lifted_product[(elements[1], elements[2])] = 0
                                else:
                                    lifted_product[(elements[1], elements[2])] = self.__zigzags_h01((bidegrees[1][0] + bidegrees[2][0], bidegrees[1][1] + bidegrees[2][1]), (elements[1]*elements[2]).basis_coefficients())
                            result = lifted_product[(elements[0], elements[1])]*elements[2] - elements[0]*lifted_product[(elements[1], elements[2])]
                            if result != 0:
                                projected_result = self.__zigzags_p(result_bidegree, result.basis_coefficients())
                                if projected_result != 0:                            
                                    dictionary[tuple(elements)] = projected_result

                if show_progress == True:
                    n += 1
                    progress = int(100*n/n_total_bidegrees)
                    if progress > previous_progress:
                        print("Progress: " + str(progress) + "%")
                    previous_progress = progress

        elif arity == 3 and operation_bidegree == (-1,-1):
            lifted_product = {}
            for bidegrees in total_bidegrees:
                if compute_symmetries or bidegrees[0][0] + bidegrees[1][0] + bidegrees[2][0] <= bidegrees[0][1] + bidegrees[1][1] + bidegrees[2][1]:
                    for elements in cartesian_product([self.zigzags_basis(bidegree) for bidegree in bidegrees]):
                        result_bidegree = (bidegrees[0][0]+bidegrees[1][0]+bidegrees[2][0]-1, bidegrees[0][1]+bidegrees[1][1]+bidegrees[2][1]-1)
                        if result_bidegree in self.bidegrees() and (compute_symmetries or (elements[2], elements[1], elements[0]) not in dictionary):
                            if (elements[0], elements[1]) not in lifted_product:
                                prod = elements[0]*elements[1]
                                if prod == 0 or (bidegrees[0][0]+bidegrees[1][0], bidegrees[0][1]+bidegrees[1][1]) not in self.bidegrees():
                                    lifted_product[(elements[0], elements[1])] = 0
                                else:
                                    lifted_product[(elements[0], elements[1])] = self.__zigzags_h11((bidegrees[0][0] + bidegrees[1][0], bidegrees[0][1] + bidegrees[1][1]), (elements[0]*elements[1]).basis_coefficients())
                            if (elements[1], elements[2]) not in lifted_product:
                                prod = elements[1]*elements[2]
                                if prod == 0 or (bidegrees[1][0]+bidegrees[2][0], bidegrees[1][1]+bidegrees[2][1]) not in self.bidegrees():
                                    lifted_product[(elements[1], elements[2])] = 0
                                else:
                                    lifted_product[(elements[1], elements[2])] = self.__zigzags_h11((bidegrees[1][0] + bidegrees[2][0], bidegrees[1][1] + bidegrees[2][1]), (elements[1]*elements[2]).basis_coefficients())
                            result = lifted_product[(elements[0], elements[1])]*elements[2] - elements[0]*lifted_product[(elements[1], elements[2])]
                            if result != 0:
                                projected_result = self.__zigzags_p(result_bidegree, result.basis_coefficients())
                                if projected_result != 0:                            
                                    dictionary[tuple(elements)] = projected_result

                if show_progress == True:
                    n += 1
                    progress = int(100*n/n_total_bidegrees)
                    if progress > previous_progress:
                        print("Progress: " + str(progress) + "%")
                    previous_progress = progress

        return dictionary

# BidifferentialBigradedCommutativeAlgebraMap
class BidifferentialBigradedCommutativeAlgebraMap(BigradedComplexMap):
    def __init__(self, domain, codomain, mapp):
        BigradedComplexMap.__init__(self, domain, codomain, mapp)

    @staticmethod
    def from_image_generators(domain, codomain, image_generators):
        domain_algebra = domain.algebra()
        codomain_algebra = codomain.algebra()
        f=GCAlgebra.Hom(domain_algebra, codomain_algebra)(image_generators)
        mapp = {bidegree:
                    Matrix(domain_algebra.base(),
                    domain.dimension(bidegree),
                    codomain.dimension(bidegree),
                    [f(b).basis_coefficients() if f(b) != 0 else [0 for _ in range(codomain.dimension(bidegree))] for b in domain_algebra.basis(bidegree)]).transpose()
                for bidegree in domain.bidegrees()}
        return BidifferentialBigradedCommutativeAlgebraMap(domain, codomain, mapp)

class BidifferentialBigradedCommutativeAlgebraExample():
    __QQi = QuadraticField(-1, 'I')

    @staticmethod
    def KodairaThurston(acs = None, names = None):
        lie_algebra = LieAlgebra(BidifferentialBigradedCommutativeAlgebraExample.__QQi, 'X,Y,Z,W', {
            ('X','Y'): {'Z':-1}
        })
        if acs == None:
            acs = Matrix(BidifferentialBigradedCommutativeAlgebraExample.__QQi,4,[
                [0,-1,0,0],
                [1,0,0,0],
                [0,0,0,-1],
                [0,0,1,0]
            ])
        if names == None:
            names = ['a', 'b', 'abar', 'bbar']
        return BidifferentialBigradedCommutativeAlgebra.from_nilmanifold(lie_algebra, acs, names)

    @staticmethod
    def Iwasawa(acs = None, names = None):
        lie_algebra = LieAlgebra(BidifferentialBigradedCommutativeAlgebraExample.__QQi, 'p,ip,q,iq,z,iz', {
            ('p','q'): {'z':1},
            ('p', 'iq'): {'iz':1},
            ('ip','q'): {'iz':1},
            ('ip', 'iq'): {'z':-1}
        })
        if acs == None:
            acs = Matrix(BidifferentialBigradedCommutativeAlgebraExample.__QQi,6,[
                [0,1,0,0,0,0],
                [-1,0,0,0,0,0],
                [0,0,0,1,0,0],
                [0,0,-1,0,0,0],
                [0,0,0,0,0,1],
                [0,0,0,0,-1,0]
            ])
        if names == None:
            names = ['a','b','c','abar','bbar','cbar']
        return BidifferentialBigradedCommutativeAlgebra.from_nilmanifold(lie_algebra, acs, names, normalization_coefficients=[1/2,1,1,1/2,1,1])

    @staticmethod
    def ST_nilmanifold():
        Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
        basis = []
        for bidegree in Iwasawa.bidegrees():
            basis += Iwasawa.algebra().basis(bidegree)
        
        generators = [basis[11], basis[15], basis[19], basis[14], basis[12], basis[22], basis[41], basis[56], basis[7]]
        return Iwasawa.subalgebra(generators)