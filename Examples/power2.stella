language core;

extend with
  #unit-type,
  #references,
  #arithmetic-operators,
  #sequencing,
  #natural-literals;

fn helper(ref : &Nat) -> fn(Nat) -> Nat {
  return
    fn (n : Nat) {
      return
        Nat::rec(n, unit, fn(i : Nat){
          return fn(r : Unit) {
            return
              Nat::rec(*ref, unit, fn(j : Nat) {
                return fn(r2 : Unit) {
                  return ref := succ(*ref)
                }
              })
          }
        });
        (*ref)
    }
}

fn exp2(n : Nat) -> Nat {
  return helper(new(1))(n)
}

fn main(n : Nat) -> Nat {
  return exp2(n)
}
