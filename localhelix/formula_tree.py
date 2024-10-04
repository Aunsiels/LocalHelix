class Node:

    def __init__(self, value, sons=None):
        self.value = value
        self.sons = sons or []

    @staticmethod
    def read_formula(formula):
        if formula.startswith("and"):
            value = "and"
            formula = formula[3:]
        elif formula.startswith("or"):
            value = "or"
            formula = formula[2:]
        elif formula.startswith("not"):
            value = "not"
            formula = formula[3:]
        else:
            return Node(formula)
        depth = 0
        current_node = Node(value)
        current = ""
        for char in formula:
            if char == "(":
                depth += 1
                if depth > 1:
                    current += char
            elif char == ")":
                depth -= 1
                if depth != 0:
                    current += char
            elif char == ",":
                if depth == 1:
                    next_node = Node.read_formula(current)
                    current_node.sons.append(next_node)
                    current = ""
                else:
                    current += char
            else:
                current += char
        if current:
            next_node = Node.read_formula(current)
            current_node.sons.append(next_node)
        return current_node

    def __str__(self):
        res = self.value
        if self.sons:
            res += "(" + ",".join(str(son) for son in self.sons) + ")"
        return res

    def check_genoset(self, genotypes):
        if self.value == "and":
            for son in self.sons:
                if not son.check_genoset(genotypes):
                    return False
            return True
        if self.value == "or":
            for son in self.sons:
                if son.check_genoset(genotypes):
                    return True
            return False
        if self.value == "not":
            if self.sons[0].check_genoset(genotypes):
                return False
            return True
        if self.value in genotypes:
            return True


if __name__ == '__main__':
    tree = Node.read_formula("and(not(gs145),rs6625163(A;A),or(rs1160312(A;A),rs1160312(A;G),and(rs201571(T;T),"
                             "rs6036025(G;G))))")
    print(tree)
    print(tree.check_genoset(["rs6625163(A;A)", "rs1160312(A;G)"]))
    print(tree.check_genoset(["rs6625163(A;A)", "rs201571(T;T)", "rs6036025(G;G)"]))
    print(tree.check_genoset(["rs6625163(A;A)"]))
    print(tree.check_genoset(["rs6625163(A;A)", "gs145"]))

