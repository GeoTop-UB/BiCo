from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace

# Auxiliary class used in the computation of the zigzags decomposition
#   ambient_space: ambient vector space where everything happens
#   outer_subspace: affine subspace of the ambient space. Points of the _PuncturedAffineSpace will be taken from this subspace
#   inner_subspaces: list of affine subspaces of the ambient space. These subspaces will be omitted when taking a point of the _PuncturedAffineSpace
class PuncturedAffineSpace():
    def __init__(self, ambient_space, outer_subspace, inner_subspaces):
        self.__ambient_space = ambient_space
        self.__ambient_dimension = ambient_space.rank()
        self.__outer_subspace = outer_subspace
        self.__inner_subspaces = []
        self.__inner_points = []
        self.__base = ambient_space.base()
        self.__is_empty = False
        try:
            self.__outer_point = outer_subspace.point()
        except:
            self.__outer_point = None
            self.__is_empty = True

        if self.__is_empty == False:
            for inner in inner_subspaces:
                intersection = inner.intersection(self.__outer_subspace)
                try:
                    p = intersection.point()
                    if intersection not in self.__inner_subspaces:
                        self.__inner_subspaces.append(intersection)
                        self.__inner_points.append(p)
                    if intersection == self.__outer_subspace:
                        self.__is_empty = True
                except:
                    pass

    @staticmethod
    def unpunctured(ambient_space, outer_subspace):
        return PuncturedAffineSpace(ambient_space, outer_subspace, [])

    @staticmethod
    def unpunctured_from_vector_space(ambient_space, outer_subspace):
        return PuncturedAffineSpace.unpunctured(ambient_space, AffineSubspace(vector([0]*ambient_space.rank()), outer_subspace))

    @staticmethod
    def empty(ambient_space):
        return PuncturedAffineSpace.unpunctured(ambient_space, None)

    @staticmethod
    def total(ambient_space):
        return PuncturedAffineSpace.unpunctured_from_vector_space(ambient_space, ambient_space)

    @staticmethod
    def zero(base, dimension):
        vector_space = VectorSpace(base, dimension)
        subspace = AffineSubspace(vector_space.zero(), vector_space.subspace([]))
        return PuncturedAffineSpace.unpunctured_from_vector_space(vector_space, subspace)

    @staticmethod
    def kernel(matrix):
        return PuncturedAffineSpace.unpunctured_from_vector_space(VectorSpace(matrix.base(), matrix.ncols()), matrix.right_kernel())

    def ambient_space(self):
        return self.__ambient_space

    def outer_subspace(self):
        return self.__outer_subspace

    def inner_subspaces(self):
        return self.__inner_subspaces

    def base(self):
        return self.__base

    def is_empty(self):
        return self.__is_empty

    def intersection(self, other):
        if self.is_empty() or other.is_empty():
            return PuncturedAffineSpace.empty(self.ambient_space())
        else:
            outer = self.outer_subspace().intersection(other.outer_subspace())
            return PuncturedAffineSpace(self.ambient_space(), outer, self.inner_subspaces() + other.inner_subspaces())

    def remove_subspace(self, remove):
        return PuncturedAffineSpace(self.ambient_space(), self.outer_subspace(), self.inner_subspaces() + [remove])

    def preimage(self, matrix):
        ambient = VectorSpace(self.base(), matrix.ncols())
        image_matrix = PuncturedAffineSpace.unpunctured_from_vector_space(self.ambient_space(), VectorSpace(self.base(), self.__ambient_dimension).subspace(matrix.columns()))
        intersection = self.intersection(image_matrix)
        if self.is_empty():
            return PuncturedAffineSpace.empty(ambient)
        else:
            try:
                x = matrix.solve_right(intersection.get_point())
                outer = AffineSubspace(x, matrix.right_kernel() + ambient.subspace([matrix.solve_right(v) for v in intersection.outer_subspace().linear_part().basis()]))
                inners = []
                for inner in intersection.inner_subspaces():
                    try:
                        y = matrix.solve_right(self.__inner_points[intersection.inner_subspaces().index(inner)])
                        inners.append(AffineSubspace(y, matrix.right_kernel() + ambient.subspace([matrix.solve_right(v) for v in inner.linear_part().basis()])))
                    except:
                        pass
                return PuncturedAffineSpace(ambient, outer, inners)
            except:
                return PuncturedAffineSpace.empty(ambient)

    def image(self, matrix):
        ambient = VectorSpace(self.base(), matrix.nrows())
        if self.is_empty():
            return PuncturedAffineSpace.empty(ambient)
        else:
            try:
                x = matrix*self.__outer_point
                outer = AffineSubspace(x, ambient.subspace([matrix*v for v in self.outer_subspace().linear_part().basis()]))
                inners = []
                for inner in self.inner_subspaces():
                    y = matrix*self.__inner_points[self.inner_subspaces().index(inner)]
                    inners.append(AffineSubspace(y, ambient.subspace([matrix*v for v in inner.linear_part().basis()])))
                return PuncturedAffineSpace(ambient, outer, inners)
            except:
                return PuncturedAffineSpace.empty(ambient)

    def get_point(self):
        if self.is_empty():
            return None
        elif self.outer_subspace().dimension() == 0:
            return self.__outer_point
        elif self.ambient_space().rank() != self.outer_subspace().linear_part().rank():
            inclusion = Matrix(self.base(), [v for v in self.outer_subspace().linear_part().basis()]).transpose()
            rank = self.outer_subspace().linear_part().rank()
            new_ambient = VectorSpace(self.base(), rank)
            new_outer = AffineSubspace(vector([0]*rank), VectorSpace(self.base(), rank))
            new_inners = []
            for inner in self.inner_subspaces():
                subspace = new_ambient.subspace([self.outer_subspace().linear_part().coordinates(v) for v in inner.linear_part().basis()])
                new_inners.append(AffineSubspace(self.outer_subspace().linear_part().coordinates(inner.point() - self.__outer_point), subspace))
            point = PuncturedAffineSpace(new_ambient, new_outer, new_inners).get_point()
            if point == None:
                return None
            else:
                return inclusion*point  + self.__outer_point
        else:
            n_varieties = len(self.inner_subspaces())
            dimension = self.ambient_space().rank()
            if dimension >= 2:
                vectors = [vector([int(i==j) for j in range(dimension)]) for i in range(dimension-2)]
                for n in range(n_varieties+1):
                    v = vector([int(i==dimension-2) for i in range(dimension)])
                    v[dimension-1] = n
                    hyperplane = AffineSubspace(vector([0]*dimension), VectorSpace(self.base(), dimension).subspace(vectors + [v]))
                    works = True
                    for inner in self.inner_subspaces():
                        if inner == hyperplane:
                            works = False
                            break

                    if works == True:
                        new_inners = []
                        for inner in self.inner_subspaces():
                            try:
                                intersection = hyperplane.intersection(inner)
                                a = intersection.point()
                                new_inners.append(intersection)
                            except:
                                pass
                        point = PuncturedAffineSpace(self.ambient_space(), hyperplane, new_inners).get_point()
                        return point
            elif dimension == 1:
                for inner in self.inner_subspaces():
                    point = inner.point() + vector([1])
                    works = True
                    for inner2 in self.inner_subspaces():
                        if point in inner2:
                            works = False
                    if works == True:
                        return point
                return self.ambient_space().zero()

# Bigraded complex
class BigradedComplex():
    def __init__(self, base, dell, delbar, names=None, latex_names=None, CHECK=True):
        self.__base = base
        self.__dimension = {}
        self.__bidegrees = []
        self.__dell = dell
        self.__delbar = delbar
        self.__delldelbar = {}
        self.__names = names
        self.__latex_names = latex_names
        self.__CHECK = CHECK

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
        self.__squares_basis = {}
        self.__zigzags_decomposition = None

        # Record the dimensions of the bigraded components
        for (p,q) in self.__dell:
            if (p,q) not in self.__dimension:
                self.__dimension[(p,q)] = self.__dell[(p,q)].ncols()
        for (p,q) in self.__delbar:
            if (p,q) not in self.__dimension:
                self.__dimension[(p,q)] = self.__delbar[(p,q)].ncols()

        for (p,q) in self.__dimension:
            if (p,q) not in self.__dell:
                if (p+1,q) not in self.__dimension:
                    self.__dell[(p,q)] = Matrix(self.base(), 0, self.__dimension[(p,q)])
                else:
                    self.__dell[(p,q)] = Matrix(self.base(), self.__dimension[(p+1,q)], self.__dimension[(p,q)])
            if (p,q) not in self.__delbar:
                if (p,q+1) not in self.__dimension:
                    self.__delbar[(p,q)] = Matrix(self.base(), 0, self.__dimension[(p,q)])
                else:
                    self.__delbar[(p,q)] = Matrix(self.base(), self.__dimension[(p,q+1)], self.__dimension[(p,q)])

        #Check if the provided bigraded complex is well-defined (i.e. if dell, delbar form a bidifferential) (in case CHECK == True)
        if self.__CHECK:
            for (p,q) in self.bidegrees():
                if (p+1,q) in self.bidegrees() and (p+2,q) in self.bidegrees():
                    if self.dell((p+1,q))*self.dell((p,q)) != 0:
                        return BaseException("The dell differential does not square to zero")
                if (p,q+1) in self.bidegrees() and (p,q+2) in self.bidegrees():
                    if self.delbar((p,q+1))*self.delbar((p,q)) != 0:
                        return BaseException("The delbar differential does not square to zero")
                if (p,q+1) in self.bidegrees() and (p+1,q) in self.bidegrees() and (p+1,q+1) in self.bidegrees():
                    if self.delbar((p+1,q))*self.dell((p,q)) + self.dell((p,q+1))*self.delbar((p,q)) != 0:
                        return BaseException("The differentials dell and delbar do not anticommute")
                if self.dell((p,q)).ncols() != self.delbar((p,q)).ncols():
                    return BaseException("The differentials dell and delbar are not well defined")
                if (p-1,q) in self.bidegrees() and self.dell((p-1,q)).nrows() != self.dell((p,q)).ncols():
                    return BaseException("The differentials dell and delbar are not well defined")
                if (p,q-1) in self.bidegrees() and self.delbar((p,q-1)).nrows() != self.dell((p,q)).ncols():
                    return BaseException("The differentials dell and delbar are not well defined")

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
                if (p,q) in self.__delbar and (p,q+1) in self.__dell and self.__delbar[(p,q)] != Matrix(self.base(), []) and self.__dell[(p,q+1)] != Matrix(self.base(), []):
                    self.__delldelbar[(p,q)] = self.__dell[(p,q+1)] * self.__delbar[(p,q)]
                else:
                    self.__delldelbar[(p,q)] = Matrix(self.base(), self.__dimension[(p+1,q+1)], self.__dimension[(p,q)])

        self.__total_degrees = []
        self.__ordered_bidegrees = {}
        for (p,q) in self.bidegrees():
            if p+q not in self.__total_degrees:
                self.__total_degrees.append(p+q)
                self.__ordered_bidegrees[p+q] = [p]
            else:
                self.__ordered_bidegrees[p+q].append(p)
        self.__total_degrees = sorted(self.__total_degrees)
        for deg in self.__total_degrees:
            self.__ordered_bidegrees[deg] = [(p,deg-p) for p in sorted(self.__ordered_bidegrees[deg])]

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
        return list(self.__dimension.keys())

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
            first = True
            result = ""
            for i in range(len(coordinates)):
                if coordinates[i] == 1:
                    if first:
                        result = str(self.names(bidegree)[i])
                        first = False
                    else:
                        result = result + " + " + str(self.names(bidegree)[i])
                elif coordinates[i] == -1:
                    if first:
                        result = "-" + str(self.names(bidegree)[i])
                        first = False
                    else:
                        result = result + " - " + str(self.names(bidegree)[i])
                elif coordinates[i] < 0:
                    if first:
                        result = str(coordinates[i]) + "*" + str(self.names(bidegree)[i])
                        first = False
                    else:
                        result = result + " - " + str(-coordinates[i]) + "*" + str(self.names(bidegree)[i])
                elif coordinates[i] > 0:
                    if first:
                        result = str(coordinates[i]) + "*" + str(self.names(bidegree)[i])
                        first = False
                    else:
                        result = result + " + " + str(coordinates[i]) + "*" + str(self.names(bidegree)[i])
            if result == "":
                result = "0"
            return result
            #return sum(c * b for (c, b) in zip(coordinates, self.names(bidegree)))

    # Return the subcomplex together with the inclusion map
    def subcomplex(self, subspace):
        valid_bidegrees = []
        for bidegree in subspace:
            if subspace[bidegree] != []:
                valid_bidegrees.append(bidegree)
        return BigradedSubcomplex({bidegree: subspace[bidegree] for bidegree in valid_bidegrees}, self, CHECK=False)

    # Return the minimal subcomplex containing the specified bidegrees
    def subcomplex_bidegrees(self, bidegrees):
        return BigradedSubcomplex({bidegree: [vector([int(i==j) for j in range(self.dimension(bidegree))]) for i in range(self.dimension(bidegree))] for bidegree in bidegrees}, self, CHECK=False)

    def total_degrees(self):
        return self.__total_degrees

    def ordered_bidegrees(self, total_degree=None):
        if total_degree == None:
            return self.__ordered_bidegrees
        else:
            return self.__ordered_bidegrees[total_degree]

    # Attach an element to the bigraded complex
    def attach_element(self, bidegree, dell, delbar, name=None):
        dell = vector(dell)
        delbar = vector(delbar)
        (p,q) = bidegree
        if (p+1,q) in self.bidegrees() and self.dell((p+1,q))*dell != 0:
            raise BaseException("The element can not be attached: dell would not square to zero.")
        elif (p,q+1) in self.bidegrees() and self.delbar((p,q+1))*delbar != 0:
            raise BaseException("The element can not be attached: delbar would not square to zero.")

        if (p+1,q) in self.bidegrees() and (p,q+1) in self.bidegrees():
            if self.delbar((p+1,q))*dell + self.dell((p,q+1))*delbar != 0:
                raise BaseException("The element can not be attached: dell and delbar would not anticommute.")
        elif (p+1,q) in self.bidegrees():
            if self.delbar((p+1,q))*dell != 0:
                raise BaseException("The element can not be attached: dell and delbar would not anticommute.")
        elif (p,q+1) in self.bidegrees():
            if self.dell((p,q+1))*delbar != 0:
                raise BaseException("The element can not be attached: dell and delbar would not anticommute.")

        new_dell = {bideg: self.dell(bideg) for bideg in self.bidegrees()}
        new_delbar = {bideg: self.delbar(bideg) for bideg in self.bidegrees()}
        if bidegree not in self.bidegrees():
            if (p+1,q) in self.bidegrees():
                new_dell[(p,q)] = Matrix(self.base(), [dell]).transpose()
            else:
                new_dell[(p,q)] = Matrix(self.base(), 0, 1)
            if (p,q+1) in self.bidegrees():
                new_delbar[(p,q)] = Matrix(self.base(), [delbar]).transpose()
            else:
                new_delbar[(p,q)] = Matrix(self.base(), 0, 1)
        else:
            if (p+1,q) not in self.bidegrees():
                new_dell[(p,q)] = Matrix(self.base(), 0, self.dimension((p,q))+1)
            else:
                values = new_dell[(p,q)].rows()
                for i in range(len(dell)):
                    values[i] = list(values[i]) + [dell[i]]
                new_dell[(p,q)] = Matrix(self.base(), values)

            if (p,q+1) not in self.bidegrees():
                new_delbar[(p,q)] = Matrix(self.base(), 0, self.dimension((p,q))+1)
            else:
                values = new_delbar[(p,q)].rows()
                for i in range(len(delbar)):
                    values[i] = list(values[i]) + [delbar[i]]
                new_delbar[(p,q)] = Matrix(self.base(), values)

        if (p-1,q) in self.bidegrees():
            values = new_dell[(p-1,q)].columns()
            if values == []:
                values = [[0]]
            else:
                for i in range(len(values)):
                    values[i] = list(values[i]) + [0]
            new_dell[(p-1,q)] = Matrix(self.base(), values).transpose()

        if (p,q-1) in self.bidegrees():
            values = new_delbar[(p,q-1)].columns()
            if values == []:
                values = [[0]]
            else:
                for i in range(len(values)):
                    values[i] = list(values[i]) + [0]
            new_delbar[(p,q-1)] = Matrix(self.base(), values).transpose()

        if self.names() != None and name != None:
            new_names = self.names()
            if (p,q) in self.bidegrees():
                new_names[(p,q)].append(name)
            else:
                new_names[(p,q)] = [name]
        else:
            new_names = None

        return BigradedComplex(self.base(), new_dell, new_delbar, names=new_names)

    def attach_random_element(self, min_degree=None, max_degree=None, min_coefficient=-5, max_coefficient=5, name=None):
        def random_element_from_vector_space(V):
            element = V.zero()
            for v in V.basis():
                coefficient = int(random()*(max_coefficient-min_coefficient)) + min_coefficient
                element += coefficient*v
            return element

        if min_degree == None or max_degree == None:
            valid_bidegree = False
            while valid_bidegree == False:
                (p,q) = self.bidegrees()[int(random()*len(self.bidegrees()))]
                if min_degree != None and (p < min_degree or q < min_degree):
                    valid_bidegree = False
                elif max_degree != None and (p > max_degree or q > max_degree):
                    valid_bidegree = False
                else:
                    valid_bidegree = True
        else:
            p = int(random()*(max_degree-min_degree))+min_degree
            q = int(random()*(max_degree-min_degree))+min_degree

        if (p+1,q) in self.bidegrees() and (p,q+1) in self.bidegrees() and (p+1,q+1) in self.bidegrees():
            dell_of_delbar_cocycles = VectorSpace(self.base(), self.dimension((p+1,q+1))).subspace([self.dell((p,q+1))*v for v in self.delbar_cocycles_raw((p,q+1)).basis()])
            delbar_of_dell_cocycles = VectorSpace(self.base(), self.dimension((p+1,q+1))).subspace([self.delbar((p+1,q))*v for v in self.dell_cocycles_raw((p+1,q)).basis()])
            subspace_dell_delbar = dell_of_delbar_cocycles.intersection(delbar_of_dell_cocycles)
            delldelbar = random_element_from_vector_space(subspace_dell_delbar)
            lift_dell = self.delbar((p+1,q)).solve_right(delldelbar)
            lift_delbar = self.dell((p,q+1)).solve_right(-delldelbar)
            subspace_dell = AffineSubspace(lift_dell, self.delbar_cocycles_raw((p+1,q))).intersection(AffineSubspace(vector([0]*self.dimension((p+1,q))), self.dell_cocycles_raw((p+1,q))))
            subspace_delbar = AffineSubspace(lift_delbar, self.dell_cocycles_raw((p,q+1))).intersection(AffineSubspace(vector([0]*self.dimension((p,q+1))), self.delbar_cocycles_raw((p,q+1))))
            dell = subspace_dell.point() + random_element_from_vector_space(subspace_dell.linear_part())
            delbar = subspace_delbar.point() + random_element_from_vector_space(subspace_delbar.linear_part())
            return self.attach_element((p,q), dell, delbar, name=name)
        elif (p+1,q) in self.bidegrees() and (p,q+1) in self.bidegrees():
            dell = random_element_from_vector_space(self.dell_cocycles_raw((p+1,q)))
            delbar = random_element_from_vector_space(self.delbar_cocycles_raw((p,q+1)))
            return self.attach_element((p,q), dell, delbar, name=name)
        elif (p+1,q) in self.bidegrees():
            if (p+1,q+1) not in self.bidegrees():
                dell = random_element_from_vector_space(self.dell_cocycles_raw((p+1,q)))
            else:
                dell = random_element_from_vector_space(self.dell_cocycles_raw((p+1,q)).intersection(self.delbar_cocycles_raw((p+1,q))))
            return self.attach_element((p,q), dell, vector([]), name=name)
        elif (p,q+1) in self.bidegrees():
            if (p+1,q+1) not in self.bidegrees():
                delbar = random_element_from_vector_space(self.delbar_cocycles_raw((p,q+1)))
            else:
                delbar = random_element_from_vector_space(self.dell_cocycles_raw((p,q+1)).intersection(self.delbar_cocycles_raw((p,q+1))))
            return self.attach_element((p,q), vector([]), delbar, name=name)
        else:
            return self.attach_element((p,q), vector([]), vector([]), name=name)

    @staticmethod
    def zero(base):
        return BigradedComplex(base, {}, {})

    @staticmethod
    def random(base, min_degree=0, max_degree=3, min_coefficient=-5, max_coefficient=5, n_generators=60, names=None):
        if names == None:
            names = ["x" + str(i+1) for i in range(n_generators)]

        partial_bicos = [BigradedComplex.zero(base)]

        for i in range(n_generators):
            partial_bicos.append(partial_bicos[-1].attach_random_element(min_degree=min_degree, max_degree=max_degree, min_coefficient=min_coefficient, max_coefficient=max_coefficient, name=names[i]))

        return partial_bicos[-1]

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
                    self.__dell_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.dell((p-1,q)).transpose())
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
                    self.__delbar_coboundaries[bidegree] = VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.delbar((p,q-1)).transpose())
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
            return proje35ction * change_basis

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
            return VectorSpace(self.base(), self.delbar_cohomology_basis(bidegree))

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

############ Zigzags ###############

    def __find_zigzag(self, already_computed={}):
        already_computed_subspaces = {}
        for bidegree in self.bidegrees():
            if bidegree in already_computed:
                already_computed_subspaces[bidegree] = AffineSubspace(vector([0]*self.dimension(bidegree)), VectorSpace(self.base(), self.dimension(bidegree)).subspace(already_computed[bidegree]))
            else:
                already_computed_subspaces[bidegree] = AffineSubspace(vector([0]*self.dimension(bidegree)), VectorSpace(self.base(), self.dimension(bidegree)).subspace([]))

        bidegrees = []
        for k in self.total_degrees():
            if bidegrees == []:
                for (p,q) in self.ordered_bidegrees(k):
                    space = PuncturedAffineSpace(VectorSpace(self.base(), self.dimension((p,q))), AffineSubspace(vector([0]*self.dimension((p,q))), self.delldelbar_cocycles_raw((p,q))), [AffineSubspace(vector([0]*self.dimension((p,q))), self.dell_coboundaries_raw((p,q)) + self.delbar_coboundaries_raw((p,q)) + already_computed_subspaces[(p,q)].linear_part())])
                    if space.is_empty() == False:
                        if bidegrees == []:
                            bidegrees.append((p,q))
                        elif (p-1,q+1) in bidegrees:
                            bidegrees.append((p,q))

        if bidegrees == []:
            return None

        (p0, q0) = bidegrees[0]

        # in this for loop we increase the length of the zigzag one unit each step
        # (the bidegree (p,q) is the other endpoint of the zigzag)
        last_cycle = False
        last = {}
        subspaces = {}
        last_p = p0-1
        zigzag = {}
        carried_intersections = {}
        if (p0,q0+1) in self.bidegrees():
            reachable_delbar = PuncturedAffineSpace(VectorSpace(self.base(), self.dimension((p0,q0+1))), AffineSubspace(vector([0]*self.dimension((p0,q0+1))), self.delbar_coboundaries_raw((p0,q0+1)).intersection(self.dell_cocycles_raw((p0,q0+1)))), [AffineSubspace(vector([0]*self.dimension((p0,q0+1))), self.delldelbar_coboundaries_raw((p0,q0+1)) + already_computed_subspaces[(p0,q0+1)].linear_part())])
            carry = reachable_delbar.preimage(self.delbar((p0,q0))).remove_subspace(already_computed_subspaces[(p0,q0)])
            if carry.is_empty():
                carry = PuncturedAffineSpace.unpunctured_from_vector_space(VectorSpace(self.base(), self.dimension((p0,q0))), self.delbar_cocycles_raw((p0,q0)))
        else:
            carry = PuncturedAffineSpace.total(VectorSpace(self.base(), self.dimension((p0,q0)))).remove_subspace(already_computed_subspaces[(p0,q0)])
        carried_intersections[(p0,q0)] = carry

        for (p,q) in bidegrees:     
            if last_p == p0-1:
                last = {}
                subspaces[(p,q)] = {(p,q): PuncturedAffineSpace(VectorSpace(self.base(), self.dimension((p,q))),
                                    AffineSubspace(vector([0]*self.dimension((p,q))), self.delldelbar_cocycles_raw((p,q))),
                                    [AffineSubspace(vector([0]*self.dimension((p,q))), self.dell_coboundaries_raw((p,q)) + self.delbar_coboundaries_raw((p,q)) + already_computed_subspaces[(p,q)].linear_part())])}
                i = p-p0

                for j in range(i):
                    subspaces[(p,q)][(p-j-1,q+j+1)] = subspaces[(p,q)][(p-j,q+j)].image(self.delbar((p-j,q+j))).preimage(self.dell((p-1-j,q+1+j))).remove_subspace(already_computed_subspaces[(p-1-j,q+1+j)])
                #carry = carry.intersection(subspaces[(p,q)][(p0,q0)])
                if (p-1,q+1) in subspaces:
                    carry = carry.intersection(subspaces[(p-1,q+1)][p0,q0])
                carried_intersections[(p,q)] = carry

                # First case
                last[(p,q)] = subspaces[(p,q)][(p,q)].intersection(PuncturedAffineSpace.unpunctured_from_vector_space(VectorSpace(self.base(), self.dimension((p,q))), self.dell_cocycles_raw((p,q)))).remove_subspace(already_computed_subspaces[(p,q)])
                for j in range(i):
                    last[(p-j-1,q+j+1)] = last[(p-j,q+j)].image(self.delbar((p-j,q+j))).preimage(self.dell((p-1-j,q+1+j))).remove_subspace(already_computed_subspaces[(p-1-j,q+1+j)])
                # i don't understand the following block of code:
                # try:
                #     a = last[(p,q)].image(self.delbar((p,q)))
                #     b = a.preimage(self.dell((p-1,q+1)))
                # except:
                #     pass
                x = (last[(p0,q0)].intersection(carry)).get_point()
                if x != None:
                    last_cycle = True
                    last_p = p
                    zigzag[(p0,q0)] = x

        # Second case
        if last_p == p0-1:
            last = {}
            (p1,q1) = bidegrees[-1]
            for i in range(p1-p0+1):
                if last_p == p0-1:
                    p = p1-i
                    q = q1+i
                    if (p+1,q) in self.bidegrees():
                        avoid = already_computed_subspaces[(p,q)].linear_part() + self.dell_cocycles_raw((p,q)) + VectorSpace(self.base(), self.dimension((p,q))).subspace([self.dell((p,q)).solve_right(v) for v in (self.delbar_coboundaries_raw((p+1,q)).intersection(self.dell_coboundaries_raw((p+1,q)))).basis()])
                        last[(p,q)] = PuncturedAffineSpace(VectorSpace(self.base(), self.dimension((p,q))), AffineSubspace(vector([0]*self.dimension((p,q))), self.delldelbar_cocycles_raw((p,q))), [AffineSubspace(vector([0]*self.dimension((p,q))), self.dell_coboundaries_raw((p,q)) + self.delbar_coboundaries_raw((p,q))), AffineSubspace(vector([0]*self.dimension((p,q))), avoid)])
                    else:
                        last[(p,q)] = PuncturedAffineSpace(VectorSpace(self.base(), self.dimension((p,q))), AffineSubspace(vector([0]*self.dimension((p,q))), self.delldelbar_cocycles_raw((p,q))), [AffineSubspace(vector([0]*self.dimension((p,q))), self.dell_coboundaries_raw((p,q)) + self.delbar_coboundaries_raw((p,q))), already_computed_subspaces[(p,q)]])
                    for j in range(p1-p0-i):
                        last[(p-j-1,q+j+1)] = last[(p-j,q+j)].image(self.delbar((p-j,q+j))).preimage(self.dell((p-j-1,q+j+1)))
                    intersection = last[(p0,q0)].intersection(carried_intersections[(p,q)])
                    x = last[(p0,q0)].intersection(carried_intersections[(p,q)]).get_point()
                    if x != None:
                        last_p = p
                        zigzag[(p0,q0)] = x


        if last_p == p0-1:
            return None

        if (p0,q0+1) in self.bidegrees():
            delbarx = self.delbar((p0,q0))*zigzag[(p0,q0)]
            if delbarx != 0:
                zigzag[(p0,q0+1)] = self.delbar((p0,q0))*zigzag[(p0,q0)]

        if (p0+1,q0) in self.bidegrees():
            dellx = self.dell((p0,q0))*zigzag[(p0,q0)]
            if dellx != 0:
                zigzag[(p0+1,q0)] = self.dell((p0,q0))*zigzag[(p0,q0)]

        for p in range(p0+1,last_p+1):
            q = p0 + q0 - p

            if p != p0:
                preimage = self.delbar((p,q)).solve_right(zigzag[(p,q+1)])
                preimage_affine = AffineSubspace(preimage, self.delbar_cocycles_raw((p,q)))
                carry = PuncturedAffineSpace.unpunctured(VectorSpace(self.base(), self.dimension((p,q))), preimage_affine)
            else:
                carry = PuncturedAffineSpace.total(VectorSpace(self.base(), self.dimension((p,q))))

            for r in range(p, last_p):
                s = p+q-r
                carry = carry.intersection(subspaces[(r,s)][(p,q)])

            zigzag[(p,q)] = carry.intersection(last[(p,q)]).get_point()

            if last_cycle == False or p != last_p:
                zigzag[(p+1,q)] = self.dell((p,q))*zigzag[(p,q)]

        return zigzag

    # Compute the zigzags decomposition
    def zigzags_decomposition(self, raw=False):
        if raw == True:
            if self.__zigzags_decomposition == None:
                try:
                    zigzags = []
                    computed = {}
                    finished = False
                    while finished == False:
                        zigzag = self.__find_zigzag(already_computed=computed)
                        if zigzag == None or zigzag == {}:
                            finished = True
                        else:
                            zigzags.append(zigzag)
                            for bidegree in zigzag:
                                if bidegree in computed:
                                    computed[bidegree].append(zigzag[bidegree])
                                else:
                                    computed[bidegree] = [zigzag[bidegree]]
                    if self._is_zigzag_decomposition(zigzags):
                        self.__zigzags_decomposition = zigzags
                        return zigzags
                    else:
                        raise BaseException("The zigzags decomposition was not computed successfully")
                except:
                    raise BaseException("The zigzags decomposition was not computed successfully")
            else:
                return self.__zigzags_decomposition
        else:
            zigzags_raw = self.zigzags_decomposition(raw=True)
            return [{bidegree: self.element(bidegree, zigzag[bidegree]) for bidegree in zigzag} for zigzag in zigzags_raw]

    # Check if a zigzag is actually a zigzag
    # Doesn't work properly... it should check if it's a MAXIMAL zigzag
    def _is_zigzag(self, zigzag):
        for (p,q) in zigzag:
            if (p+1,q) in zigzag:
                if self.dell((p,q))*zigzag[(p,q)] != zigzag[(p+1,q)]:
                    return False
            else:
                if self.dell((p,q))*zigzag[(p,q)] != 0:
                    return False
            if (p,q+1) in zigzag:
                if self.delbar((p,q))*zigzag[(p,q)] != zigzag[(p,q+1)]:
                    return False
            else:
                if self.delbar((p,q))*zigzag[(p,q)] != 0:
                    return False
        return True

    # Check if the given zigzag decomposition forms a zigzag decomposition
    def _is_zigzag_decomposition(self, zigzags):
        subcomplex = self.subcomplex({})
        for zigzag in zigzags:
            zigzag_subcomplex = self.subcomplex({bidegree: [zigzag[bidegree]] for bidegree in zigzag})
            if zigzag_subcomplex.is_zigzag() == False:
                return False
            elif sum(subcomplex.intersection(zigzag_subcomplex).dimension(bidegree) for bidegree in self.bidegrees()) != 0:
                return False
            else:
                subcomplex = zigzag_subcomplex.sum(subcomplex)
        for bidegree in self.bidegrees():
            if subcomplex.dimension(bidegree) != self.reduced_aeppli_cohomology_raw(bidegree).rank() + self.bottchern_cohomology_raw(bidegree).rank():
                return False
        return True

    # Return the (a) subcomplex of zigzags
    def zigzags_subcomplex(self):
        basis = {bidegree: [] for bidegree in self.bidegrees()}
        for zigzag in self.zigzags_decomposition(raw=True):
            for bidegree in zigzag:
                basis[bidegree].append(zigzag[bidegree])
        return self.subcomplex(basis)

    # Return the (a) basis of zigzags
    def zigzags_basis(self, bidegree=None, raw=False):
        if raw == True:
            zigzags_subcomplex = self.zigzags_subcomplex()
            if bidegree == None:
                return {bidegree: [v for v in zigzags_subcomplex.subspace(bidegree).basis()] for bidegree in self.bidegrees()}
            else:
                return [v for v in zigzags_subcomplex.subspace(bidegree).basis()]
        else:
            basis = self.zigzags_basis(bidegree=bidegree, raw=True)
            if bidegree == None:
                return {bidegree: [self.element(bidegree, v) for v in basis[bidegree]] for bidegre in self.bidegrees()}
            else:
                return [self.element(bidegree, v) for v in basis]


    # Compute the shapes of the zigzags in the zigzag decomposition
    # The shapes are expressed in the notation from "On the structure of double complexes", section 2 (by Jonas Stelzig)
    def zigzags_shapes(self, raw=False):
        zigzags = self.zigzags_decomposition(raw=raw)
        zigzags_even = {}
        zigzags_odd = {}
        for zigzag in zigzags:
            total_degrees = {}
            length = len(zigzag)
            top_corner_p = None
            top_corner_q = None
            bot_corner_p = None
            bot_corner_q = None
            for (p,q) in zigzag:
                if p+q in total_degrees:
                    total_degrees[p+q] += 1
                else:
                    total_degrees[p+q] = 1
                if top_corner_p == None or p < top_corner_p or q > top_corner_q:
                    top_corner_p = p
                    top_corner_q = q
                if bot_corner_p == None or p > bot_corner_p or q < bot_corner_q:
                    bot_corner_p = p
                    bot_corner_q = q

            # Even zigzag
            if length % 2 == 0:
                if (top_corner_p + 1, top_corner_q) in zigzag:
                    i = 1
                else:
                    i = 2
                r = length/2
                if (top_corner_p, top_corner_q, i, length) in zigzags_even:
                    zigzags_even[(top_corner_p, top_corner_q, i, length)].append(zigzag[(top_corner_p, top_corner_q)])
                else:
                    zigzags_even[(top_corner_p, top_corner_q, i, length)] = [zigzag[(top_corner_p, top_corner_q)]]
            # Odd zigzag
            else:
                # Dot
                if len(total_degrees) == 1:
                    if (top_corner_p, top_corner_q, top_corner_p+top_corner_q) in zigzags_odd:
                        zigzags_odd[(top_corner_p, top_corner_q, top_corner_p+top_corner_q)].append(zigzag[(top_corner_p, top_corner_q)])
                    else:
                        zigzags_odd[(top_corner_p, top_corner_q, top_corner_p+top_corner_q)] = [zigzag[(top_corner_p, top_corner_q)]]
                # Proper zigzag
                else:
                    [deg1, deg2] = list(total_degrees.keys())
                    if total_degrees[deg1] > total_degrees[deg2]:
                        d = deg1
                        if deg2 > deg1:
                            data = (top_corner_p, bot_corner_q, d)
                            
                        else:
                            data = (bot_corner_p, top_corner_q, d)
                    else:
                        d = deg2
                        if deg1 > deg2:
                            data = (top_corner_p, bot_corner_q, d)
                        else:
                            data = (bot_corner_p, top_corner_q, d)
                    if data in zigzags_odd:
                        zigzags_odd[data].append(zigzag[(top_corner_p, top_corner_q)])
                    else:
                        zigzags_odd[data] = [zigzag[(top_corner_p, top_corner_q)]]

        return (zigzags_even, zigzags_odd)

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

######################## Squares ###########################

    # Squares
    def squares_basis(self, bidegree=None, raw=False):
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
        if bidegree != None:
            if bidegree not in self.__squares_basis:
                zigzags_basis = self.zigzags_basis(bidegree, raw=True)
                self.__squares_basis[bidegree] = []
                for i in range(self.dimension(bidegree)):
                    v = vector([int(i==j) for j in range(self.dimension(bidegree))])
                    if v not in VectorSpace(self.base(), self.dimension(bidegree)).subspace(zigzags_basis + self.__squares_basis[bidegree]):
                        self.__squares_basis[bidegree] += [v]
            if raw == True or self.names() == None:
                return self.__squares_basis[bidegree]
            else:
                return [self.element(bidegree, b) for b in self.squares_basis(bidegree, raw=True)]
        else:
            for bidegree in self.bidegrees():
                self.squares_basis(bidegree=bidegree, raw=True)
            if raw == True:
                return self.__squares_basis
            else:
                return {bidegree: [self.element(bidegree, b) for b in self.__squares_basis[bidegree]] for bidegree in self.bidegrees()}

    def squares(self, bidegree=None, raw=False):
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
        if bidegree != None:
            if raw == True or self.__names == None:
                return VectorSpace(self.base(), self.dimension(bidegree)).subspace(self.squares_basis(bidegree, raw=True))
            else:
                return VectorSpace(self.base(), self.squares_basis(bidegree))
        else:
            return {bidegree: self.squares(bidegree=bidegree, raw=raw) for bidegree in self.bidegrees()}

    # Squares decomposition of the bigraded complex
    # Squares contain: vertex, dell(vertex), delbar(vertex), dell delabr(vertex)
    def squares_decomposition(self, raw=False):
        if raw == True:
            squares = self.squares_basis(raw=True)
            decomposition = []
            computed_squares = {bidegree: [] for bidegree in squares}
            zigzags = self.zigzags_basis(raw=True)
            for (p,q) in squares:
                for vertex in squares[(p,q)]:
                    if self.delldelbar((p,q))*vertex != 0 and vertex not in VectorSpace(self.base(), self.dimension((p,q))).subspace(zigzags[(p,q)] + computed_squares[(p,q)]):
                        dell = self.dell((p,q))*vertex
                        delbar = self.delbar((p,q))*vertex
                        delldelbar = self.delldelbar((p,q))*vertex
                        square = {(p,q): vertex, (p+1,q): dell, (p,q+1): delbar, (p+1,q+1): delldelbar}
                        computed_squares[(p,q)].append(vertex)
                        computed_squares[(p+1,q)].append(dell)
                        computed_squares[(p,q+1)].append(delbar)
                        computed_squares[(p+1,q+1)].append(delldelbar)
                        decomposition.append(square)
            return decomposition                    
        else:
            return [{bidegree: self.element(bidegree, square[bidegree]) for bidegree in square} for square in self.squares_decomposition(raw=True)]

############# Frölicher spectral sequence #############

    def total_degree_to_bidegree(self, total_degree, element):
        result = {}
        index = 0
        for bidegree in self.ordered_bidegrees(total_degree):
            result[bidegree] = vector(element[index : index + self.dimension(bidegree)])
            index = index + self.dimension(bidegree)
        return result

    def bidegree_to_total_degree(self, element):
        if element == {}:
            return 0
        else:
            (p,q) = list(element.keys())[0]
            total_degree = p+q
            total_element = []
            for bidegree in self.ordered_bidegrees(total_degree):
                if bidegree in element:
                    total_element = total_element + list(element[bidegree])
                else:
                    total_element = total_element + [0]*self.dimension(bidegree)
            return vector(total_element)

    def total_dimension(self, total_degree=None):
        if total_degree != None:
            return sum(self.dimension(bidegree) for bidegree in self.ordered_bidegrees(total_degree))
        else:
            return {total_degree: self.total_dimension(total_degree=total_degree) for total_degree in self.total_degrees()}

    def total_dell(self, total_degree):
        matrix = []
        for (p,q) in self.ordered_bidegrees(total_degree):
            for column in self.dell((p,q)).columns():
                matrix.append(self.bidegree_to_total_degree({(p+1,q): column}))
        return Matrix(self.base(), matrix).transpose()

    def total_delbar(self, total_degree):
        matrix = []
        for (p,q) in self.ordered_bidegrees(total_degree):
            for column in self.delbar((p,q)).columns():
                matrix.append(self.bidegree_to_total_degree({(p,q+1): column}))
        return Matrix(self.base(), matrix).transpose()

    def total_differential(self, total_degree):
        if total_degree in self.total_degrees():
            if total_degree+1 in self.total_degrees():
                return (self.total_dell(total_degree) + self.total_delbar(total_degree))
            else:
                return Matrix(self.base(), 0, self.total_dimension(total_degree))
        else:
            if total_degree+1 in self.total_degrees():
                return Matrix(self.base(), self.total_dimension(total_degree+1), 0)
            else:
                return Matrix(self.base(), 0)
        

    def total_coboundaries(self, total_degree):
        if total_degree not in self.total_degrees():
            return VectorSpace(self.base(), 0)
        elif total_degree-1 not in self.total_degrees():
            return VectorSpace(self.base(), self.total_dimension(total_degree)).subspace([])
        else:
            return VectorSpace(self.base(), self.total_dimension(total_degree)).subspace(self.total_differential(total_degree-1).columns())

    def total_cocycles(self, total_degree):
        if total_degree not in self.total_degrees():
            return VectorSpace(self.base(), 0)
        else:
            return self.total_differential(total_degree).right_kernel()

    def horizontal_filtration(self, stage, total_degree):
        if total_degree not in self.total_degrees():
            return VectorSpace(self.base(), 0)
        else:
            basis = []
            for (p,q) in self.ordered_bidegrees(total_degree):
                if p >= stage:
                    for i in range(self.dimension((p,q))):
                        v = vector([int(i==j) for j in range(self.dimension((p,q)))])
                        basis.append(self.bidegree_to_total_degree({(p,q): v}))
            return VectorSpace(self.base(), self.total_dimension(total_degree)).subspace(basis)

    def spectral_sequence_cocycles(self, r, p, total_degree):
        if total_degree in self.total_degrees():
            if total_degree+1 not in self.total_degrees():
                return self.horizontal_filtration(p, total_degree)
            else:
                intersection = VectorSpace(self.base(), self.total_dimension(total_degree+1)).subspace([self.total_differential(total_degree)*v for v in self.horizontal_filtration(p, total_degree).basis()]).intersection(self.horizontal_filtration(p+r, total_degree+1))
                return self.total_cocycles(total_degree).intersection(self.horizontal_filtration(p, total_degree)) + VectorSpace(self.base(), self.total_dimension(total_degree)).subspace(self.total_differential(total_degree).solve_right(v) for v in intersection.basis())
        else:
            return VectorSpace(self.base(), 0)

    def spectral_sequence_coboundaries(self, r, p, total_degree):
        if total_degree in self.total_degrees():
            return self.horizontal_filtration(p, total_degree).intersection(VectorSpace(self.base(), self.total_dimension(total_degree)).subspace([self.total_differential(total_degree-1)*v for v in self.horizontal_filtration(p-r, total_degree-1).basis()]))
        else:
            return VectorSpace(self.base(), 0)

    def spectral_sequence_basis(self, r, p, q, raw=False):
        if raw == True or self.names() == None:
            return (self.spectral_sequence_cocycles(r, p, p+q)/(self.spectral_sequence_cocycles(r-1,p+1,p+q)+self.spectral_sequence_coboundaries(r-1,p,p+q))).basis()
        else:
            spectral_sequence_raw = self.spectral_sequence(r, p, q, raw=True)
            lifted_basis = [self.total_degree_to_bidegree(p+q, spectral_sequence_raw.lift(b)) for b in spectral_sequence_raw.basis()]
            basis = []
            for b in lifted_basis:
                element = ""
                for bidegree in b:
                    if b[bidegree] != 0:
                        element = element + str(self.element(bidegree, b[bidegree])) + " + "
                basis.append('[{}]'.format(element[:-3]))
            return basis

    def spectral_sequence(self, r, p, q, raw=False):
        if raw == True or self.names() == None:
            return self.spectral_sequence_cocycles(r, p, p+q)/(self.spectral_sequence_cocycles(r-1,p+1,p+q)+self.spectral_sequence_coboundaries(r-1,p,p+q))
        else:
            return VectorSpace(self.base(), self.spectral_sequence_basis(r,p,q,raw=False))

############# Ascii art ################

    # Method used to produce the asii art tables for the cohomologies
    # Data is a dictionary where each (p,q)-entry contains a basis for the (p,q)-th component
    # n_generators_row is an integer parameter that determines the maximum number of generators that will appear per row on each cell of the table. If it is set to -1, there will be no limit
    def __ascii_art_table(self, data, n_generators_row):
        width = {}
        height = {}
        total_width = 0
        string = {}

        for bidegree in self.bidegrees():
            if bidegree not in data:
                data[bidegree] = []

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

    def _ascii_art_zigzags_decomposition(self, raw=False):
        zigzags = self.zigzags_decomposition(raw=raw)
        count = 1
        for zigzag in zigzags:
            print("Zigzag #" + str(count) + ": \n")
            self.__ascii_art_table({bidegree: [zigzag[bidegree]] for bidegree in zigzag}, n_generators_row=-1)
            print("\n\n")
            count += 1

    def _ascii_art_dots(self, n_generators_row=-1):
        zigzags = self.zigzags_decomposition()
        dots = {}
        for zigzag in zigzags:
            if len(zigzag) == 1:
                for bidegree in zigzag:
                    if bidegree in dots:
                        dots[bidegree].append(zigzag[bidegree])
                    else:
                        dots[bidegree] = [zigzag[bidegree]]
        self.__ascii_art_table(dots, n_generators_row=n_generators_row)

    def _ascii_art_spectral_sequence(self, r, n_generators_row=-1):
        self.__ascii_art_table({bidegree: self.spectral_sequence_basis(r, bidegree[0], bidegree[1]) for bidegree in self.bidegrees()}, n_generators_row=n_generators_row)

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

# BigradedSubcomplex
class BigradedSubcomplex(BigradedComplex):
    def __init__(self, basis, parent, CHECK=False):
        self.__parent = parent

        new_basis = {bidegree: [] for bidegree in basis}
        for (p,q) in basis:
            for x in basis[(p,q)]:
                new_basis[(p,q)].append(x)
                dellx = self.__parent.dell((p,q))*x
                delbarx = self.__parent.delbar((p,q))*x
                delldelbarx = self.__parent.delldelbar((p,q))*x
                if dellx != 0:
                    if (p+1,q) not in new_basis:
                        new_basis[(p+1,q)] = [dellx]
                    else:
                        new_basis[(p+1,q)].append(dellx)
                if delbarx != 0:
                    if (p,q+1) not in new_basis:
                        new_basis[(p,q+1)] = [delbarx]
                    else:
                        new_basis[(p,q+1)].append(delbarx)
                if delldelbarx != 0:
                    if (p+1,q+1) not in new_basis:
                        new_basis[(p+1,q+1)] = [delldelbarx]
                    else:
                        new_basis[(p+1,q+1)].append(delldelbarx)

        self.__subspace = {bidegree: VectorSpace(self.__parent.base(), self.__parent.dimension(bidegree)).subspace(new_basis[bidegree]) for bidegree in new_basis}
        dell = {}
        delbar = {}
        for (p,q) in self.__subspace:
            if (p+1,q) not in self.__subspace:
                dell[(p,q)] = Matrix(self.__parent.base(), 0, self.__subspace[(p,q)].rank())
            else:
                dell[(p,q)] = Matrix(self.__parent.base(), [self.__subspace[(p+1,q)].coordinates(self.__parent.dell((p,q))*x) for x in self.__subspace[(p,q)].basis()]).transpose()
            if (p,q+1) not in self.__subspace:
                delbar[(p,q)] = Matrix(self.__parent.base(), 0, self.__subspace[(p,q)].rank())
            else:
                delbar[(p,q)] = Matrix(self.__parent.base(), [self.__subspace[(p,q+1)].coordinates(self.__parent.delbar((p,q))*x) for x in self.__subspace[(p,q)].basis()]).transpose()

        if self.__parent.names() != None:
            names = {bidegree: ["(" + str(self.__parent.element(bidegree, x)) + ")" for x in self.__subspace[bidegree].basis()] for bidegree in self.__subspace}
        else:
            names = None

        BigradedComplex.__init__(self, self.__parent.base(), dell, delbar, names=names, CHECK=CHECK)

        inclusion_matrix = {bidegree: Matrix(self.__parent.base(), [x for x in self.__subspace[bidegree].basis()]).transpose() for bidegree in self.__subspace}
        self.__inclusion = BigradedComplexMap(self, self.__parent, inclusion_matrix)

    def subspace(self, bidegree=None):
        if bidegree == None:
            return self.__subspace
        else:
            if bidegree in self.bidegrees():
                return self.__subspace[bidegree]
            else:
                return VectorSpace(self.base(), self.parent().dimension(bidegree)).subspace([])

    def inclusion(self):
        return self.__inclusion

    def parent(self):
        return self.__parent

    def parent_element(self, bidegree, element, raw=True):
        return self.__inclusion(bidegree, element, raw=raw)

    def element(self, bidegree, element):
        return self.parent().element(bidegree, self.parent_element(bidegree, element))

    def intersection(self, subcomplex):
        intersection_subspace = {}
        for bidegree in self.bidegrees():
            if bidegree in subcomplex.bidegrees():
                intersection_subspace[bidegree] = self.subspace(bidegree=bidegree).intersection(subcomplex.subspace(bidegree=bidegree))
        return self.parent().subcomplex({bidegree: intersection_subspace[bidegree].basis() for bidegree in intersection_subspace})

    def sum(self, subcomplex):
        sum_subspace = {}
        for bidegree in self.bidegrees():
            if bidegree in subcomplex.bidegrees():
                sum_subspace[bidegree] = self.subspace(bidegree=bidegree) + subcomplex.subspace(bidegree=bidegree)
            else:
                sum_subspace[bidegree] = self.subspace(bidegree=bidegree)
        for bidegree in subcomplex.bidegrees():
            if bidegree not in sum_subspace:
                if bidegree in self.bidegrees():
                    sum_subspace[bidegree] = self.subspace(bidegree=bidegree) + subcomplex.subspace(bidegree=bidegree)
                else:
                    sum_subspace[bidegree] = subcomplex.subspace(bidegree=bidegree)
        return self.parent().subcomplex({bidegree: sum_subspace[bidegree].basis() for bidegree in sum_subspace})

    #Check if the subcomplex is a zigzag
    def is_zigzag(self):
        # TODO: remains to be checked if the zigzag is part of a square!
        # (isn't it already done?)
        if len(self.bidegrees()) == 0:
            return False
        elif len(self.bidegrees()) == 1:
            return self.dimension(self.bidegrees()[0]) == 1
        else:
            total_degrees = self.total_degrees()
            bidegrees = self.ordered_bidegrees()
            if len(total_degrees) != 2:
                return False
            if total_degrees[1] != total_degrees[0]+1:
                return False
            for (p,q) in self.bidegrees():
                if self.dimension((p,q)) != 1:
                    return False
            (r,s) = bidegrees[total_degrees[0]][0]
            (u,v) = bidegrees[total_degrees[0]][-1]
            for (p,q) in bidegrees[total_degrees[0]]:
                if (p,q) == (r,s):
                    if (r,s) != (u,v):
                        if (r+1,s) not in self.bidegrees() or self.dell((r,s)) == 0:
                            return False
                elif (p,q) == (u,v):
                    if (u,v+1) not in self.bidegrees() or self.delbar((u,v)) == 0:
                        return False
                else:
                    if (p+1,q) not in self.bidegrees() or (p,q+1) not in self.bidegrees() or self.dell((p,q)) == 0 or self.delbar((p,q)) == 0:
                        return False
            for (p,q) in bidegrees[total_degrees[0]]:
                if (self.parent().dell_coboundaries_raw((p,q)) + self.parent().delbar_coboundaries_raw((p,q))).intersection(self.subspace((p,q))).rank() != 0:
                    return False
            return True

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

        zero_dell_differential = True
        for x in dell_dictionary:
            if dell_dictionary[x] != 0:
                zero_dell_differential = False
                break
        zero_delbar_differential = True
        for x in delbar_dictionary:
            if delbar_dictionary[x] != 0:
                zero_delbar_differential = False
                break

        if zero_dell_differential == False:
            dell_differential = self.__algebra.differential(dell_dictionary)
        if zero_delbar_differential == False:
            delbar_differential = self.__algebra.differential(delbar_dictionary)
        dell_matrix = {}
        delbar_matrix = {}
        bigraded_basis = {}
        for p in range(self.__min_deg, self.__max_deg+1):
            for q in range(self.__min_deg, self.__max_deg+1):
                basis = self.__algebra.basis((p,q))
                if basis != []:
                    bigraded_basis[(p,q)] = basis

                    # Compute the differential matrices
                    if zero_dell_differential:
                        dell_matrix[(p,q)] = Matrix(self.__algebra.base(), len(self.__algebra.basis((p+1,q))), len(basis))
                    else:
                        dell_matrix[(p,q)] = dell_differential.differential_matrix_multigraded((p,q)).transpose()
                    if zero_delbar_differential:
                        delbar_matrix[(p,q)] = Matrix(self.__algebra.base(), len(self.__algebra.basis((p,q+1))), len(basis))
                    else:
                        delbar_matrix[(p,q)] = delbar_differential.differential_matrix_multigraded((p,q)).transpose()

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

    # TODO: does not work
    # Orbifold defined in "Dolbeault and Bott-Chern formalities: deformations and dell-delbar lemma", by Tomasso Sferruzza and Adriano Tomassini
    # @staticmethod
    # def SferruzzaTomassini_orbifold():
    #     Iwasawa = BidifferentialBigradedCommutativeAlgebraExample.Iwasawa()
    #     basis = []
    #     for bidegree in Iwasawa.bidegrees():
    #         basis += Iwasawa.algebra().basis(bidegree)
        
    #     generators = [basis[11], basis[15], basis[19], basis[14], basis[12], basis[22], basis[41], basis[56], basis[7]]
    #     return Iwasawa.subalgebra(generators)
