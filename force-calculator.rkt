#lang racket
(require "declarations.rkt")
(provide buildTree calcForces moveparticles)
(define (buildTree initialArea particles)
  (define (make-subtrees)
    (let* ([x1 (bbox-llx initialArea)]
           [y1 (bbox-lly initialArea)]
           [x2 (vec-x (center initialArea))]
           [y2 (vec-y (center initialArea))]
           [x3 (bbox-rux initialArea)]
           [y3 (bbox-ruy initialArea)]
           [ul-area (bbox x1 y2 x2 y3)]
           [ur-area (bbox x2 y2 x3 y3)]
           [ll-area (bbox x1 y1 x2 y2)]
           [lr-area (bbox x2 y1 x3 y2)]
           [particles-ul-area (lc x : x <- particles @ (check-in-square x x1 y2 x2 y3))]
           [particles-ur-area (lc x : x <- particles @ (check-in-square x x2 y2 x3 y3))]
           [particles-ll-area (lc x : x <- particles @ (check-in-square x x1 y1 x2 y2))]
           [particles-lr-area (lc x : x <- particles @ (check-in-square x x2 y1 x3 y2))])
      (filter non-empty-mass? (list (buildTree ul-area particles-ul-area)
                               (buildTree ur-area particles-ur-area)
                               (buildTree ll-area particles-ll-area)
                               (buildTree lr-area particles-lr-area)))))
  (define (check-in-square particle-a llx lly rux ruy)
    (and (and (>= (vec-x (particle-posn particle-a)) llx) (< (vec-x (particle-posn particle-a)) rux))
         (and (>= (vec-y (particle-posn particle-a)) lly) (< (vec-y (particle-posn particle-a)) ruy))))
  (cond [(null? particles) (gnode 0 (center initialArea) '())]
        [(singleton particles) (gnode (particle-mass (car particles)) (particle-posn (car particles)) '())]
        [#t (let ((list-of-subtrees (make-subtrees)))
              (define (total-mass-moment)
                (vec (foldr + 0 (map (lambda (tree) (* (vec-x (gnode-posn tree)) (gnode-mass tree))) list-of-subtrees))
                     (foldr + 0 (map (lambda (tree) (* (vec-y (gnode-posn tree)) (gnode-mass tree))) list-of-subtrees))))
              (define (total-mass)
                (foldr + 0 (map gnode-mass list-of-subtrees)))
              (define (com)
                (vec (/ (vec-x (total-mass-moment)) (total-mass))
                     (/ (vec-y (total-mass-moment)) (total-mass))))
              (gnode (total-mass) (com) list-of-subtrees))]))
(define (center area)
  (vec (/ (+ (bbox-llx area)
                            (bbox-rux area)) 2)
                      (/ (+ (bbox-lly area)
                            (bbox-ruy area)) 2)))
(define (non-empty-mass? node)
  (not (= (gnode-mass node) 0)))
(define (calcForces initialArea tree particles)
  (define (force-on-one-helper particle-1 s tree)
    (define (near? posn1 posn2)
    (<= (/ (dist posn1 posn2) s) 2))
    (cond [(and (near? (gnode-posn tree) (particle-posn particle-1)) (not (null? (gnode-subtrees tree))))
           (sum-of-forces (map (lambda (subtree) (force-on-one-helper particle-1 (/ s 2) subtree)) (gnode-subtrees tree)))]
          [#t (newton-force (particle-mass particle-1)
                            (gnode-mass tree)
                            (particle-posn particle-1)
                            (gnode-posn tree))]))
  (define (force-on-one particle-1 tree)
    (force-on-one-helper particle-1 (- (bbox-rux initialArea) (bbox-llx initialArea)) tree))
  (define (sum-of-forces list-of-forces)
    (vec (foldr + 0 (map vec-x list-of-forces))
         (foldr + 0 (map vec-y list-of-forces))))
  (define (newton-force mass1 mass2 posn1 posn2)
    (if (= (dist posn1 posn2) 0) (vec 0 0)
    (vec (/ (* g mass1 mass2 (- (vec-x posn2) (vec-x posn1))) (expt (dist posn1 posn2) 3))
         (/ (* g mass1 mass2 (- (vec-y posn2) (vec-y posn1))) (expt (dist posn1 posn2) 3)))))
  (define (dist posn1 posn2)
    (expt (+ (expt (- (vec-x posn2) (vec-x posn1)) 2) (expt (- (vec-y posn2) (vec-y posn1)) 2)) 0.5))
  (map (lambda (a) (force-on-one a tree)) particles))
(define (moveparticles particles forces)
  (zipwith movement particles forces))
(define (movement particle-a force-on-a)
  (particle (particle-mass particle-a)
            (vec (+ (vec-x (particle-posn particle-a)) (* (vec-x (particle-velocity particle-a)) timeslice) (/ (* (vec-x (accleration particle-a force-on-a)) timeslice timeslice) 2))
                 (+ (vec-y (particle-posn particle-a)) (* (vec-y (particle-velocity particle-a)) timeslice) (/ (* (vec-y (accleration particle-a force-on-a)) timeslice timeslice) 2)))
            (vec (+ (vec-x (particle-velocity particle-a)) (* (vec-x (accleration particle-a force-on-a)) timeslice))
                 (+ (vec-y (particle-velocity particle-a)) (* (vec-y (accleration particle-a force-on-a)) timeslice)))))
(define (accleration particle-1 force-on-1)
  (vec (/ (vec-x force-on-1) (particle-mass particle-1))
       (/ (vec-y force-on-1) (particle-mass particle-1))))