
(defvar test nil "variable to hold the point data")
(defvar test-metadata "centroids, varainces and angle of created clusters")
(defparameter *distribution* nil "integer uniform distribution or guassian")
(defparameter *min-variance* 5 "Minimum variance for each dimension in a cluster")
(defparameter *variance-range* 200 "Range of variance for each dimension in a cluster")
(defparameter *num-clusters* 6 "number of clusters to generate")
(defparameter *num-points* 1000 "number of points to generate")

(defun gaussian-random (mean variance)
  "Generates a random number under a gaussian distribution with the given mean and variance (using the Box-Muller-Marsaglia method)"
  (let (x y (w 0))
    (while (not (and (< 0 w) (< w 1))) '()
      (setf x (- (random 2.0) 1.0))
      (setf y (- (random 2.0) 1.0))
      (setf w (+ (* x x) (* y y))))
    (+ mean (* x (sqrt variance) (sqrt (* -2 (/ (log w) w)))))))

(defun random? (&optional (prob 0.5))
  "Tosses a coin of prob probability of coming up heads,
then returns t if it's heads, else nil."
  (< (random 1.0) prob))

(defun integer-random (mean variance)
  "Generates a random integer with uniform distribution with a range of 2x variance."
  (funcall (if (random?) #'+ #'-) mean (random variance))) 

(defun create-centroids (num-clusters)
  "Generates n-dimensional centroids to represent the mean of the clusters"
  (generate-list num-clusters
		 (lambda () (generate-list *point-dimensions* (lambda ()  (- (random (* 2 *max-value*)) *max-value*))
					   t))
		 t))

(defun random-variances ()
  "Creates a vector of random variances for generating clusters in *point-dimensions* space"
  (generate-list *point-dimensions* (lambda () (+ *min-variance* (random *variance-range*)))))

(defun rotate-about-origin (angle x y &optional (prototype 1.00000))
  "Rotates (x, y) about the origin. x = x1cos - y1sin ; y = y1cos + x1sin"
  (list (float (- (* x (cos angle)) (* y (sin angle))) prototype) (float (+ (* x (sin angle)) (* y (cos angle))) prototype)))

(defun 2d-rotation (point angle centroid)
  "if *point-dimensions* is 2, will rotate the new point angle 'degrees around the centroid"
  (if (/= *point-dimensions* 2) point ; return point if its not 2D
      (let* ((delta (mapcar #'- point centroid)) ; move to origin
	     (rotated (apply #'rotate-about-origin angle delta)) ; rotate point
	     (result (mapcar #'+ rotated centroid))) ; move back to centroid
	(display 4 "rotate -> delta=~A rotated=~A" delta rotated)
	result)))
	  
(defun create-point (cluster-data)
  "Given cluster data of the form '(centroid variances angle), creates a point near centroid with 
given set of variances (one for each dimension) rotated 'angle degrees"
  (let ((centroid (first cluster-data)) ; unpack cluster-data
	(variances (second cluster-data)) 
	(angle (third cluster-data))
	result component)
    (dotimes (i *point-dimensions* result) ; iterate over x, y, z, ...
      (setf component (funcall *distribution* (elt centroid i) (elt variances i))) ; result i = dist(xi vi)
      (setf result (append result (list component)))) ; append point to result
    (2d-rotation result angle centroid))) ; rotate new point if it is 2d

(defun random-clusters (&key (num-clusters *num-clusters*) (num-points *num-points*))
  "creates a data set with variable cluster sizes"
  (let* ((centroids (create-centroids num-clusters)) ; create centroids to base clusters around
	 (variances (generate-list num-clusters #'random-variances)) ; create a dimension-variances list for each centroid
	 (angles (generate-list num-clusters (lambda () (random pi)))) ; create an angle (up to pi) for each cluster
	 (cluster-data (mapcar #'list centroids variances angles))) ; pack data needed to make a cluster 
    (display 1 "~%Centroid, Variances, Angle = ~A~%" cluster-data)
    (setf test-metadata cluster-data)
    (generate-list num-points (lambda () (create-point (random-elt cluster-data)))))) ; pick cluster to add point to 
  
(defun test-create-centriods ()
  (dotimes (i 5)
    (display 0 "~%Dimensions=~D" (1+ i))
    n(setf *point-dimensions* (1+ i))
    (dotimes (j 5)
      (print (create-centroids (1+ j))))))

(defun test-create-point (variances)
  (let (centroid)
    (dotimes (i 5)
      (display 0 "~%Dimensions=~D" (1+ i))
      (setf *point-dimensions* (1+ i))
      (setf centroid (first (create-centroids 1)))
      (display 0 "~%Centroid  = ~A" centroid)
      (dotimes (j 5)
	(print (create-point (list centroid variances 0)))))))
    
(defun test-random-clusters ()
  (let ((tmp *variance-range*)
	(tmp2 *min-variance*)
	clusters points)
    (setf *variance-range* 1)
    (dotimes (i 3)
      (setf clusters (1+ i))
      (display 0 "~%~%~%----------------------- Num-Clusters ~D -----------------------------" clusters)
      (dotimes (j 3)
	(setf points  (* 4 (1+ j)))
	(display 0 "~%~%.......... Num-Points ~D............" points)
	(dotimes (k 3)
	  (setf *min-variance* (1+ (* 5 k)))
	  (display 0 "~%Variance ~D" *min-variance*)
	  (print (random-clusters clusters points)))))
    (setf *variance-range* tmp *min-variance* tmp2)))	

(defun test-lp-distance ()
  (let (centroid variance point1 point2 p)
    (dotimes (i 4)
      (setf *point-dimensions* (1+ i))
      (setf centroid (make-list *point-dimensions* :initial-element 0))
      (display 0 "~%~%~%----------------------- Dimensions ~D -----------------------------~%Centriod = ~A" *point-dimensions* centroid)
      (dotimes (j 3)
	(setf variance (* 4 (1+ j)))
	(setf point1 (create-point (list centroid variance 0)) point2 (create-point (list centroid variance 0)))
	(display 0 "~%~%.......... Variance ~D............~%point1=~A~%point2=~A~%" variance point1 point2)
	(dotimes (k 3)
	  (setf p (1+ k))
	  (display 0 "~%p=~D -> ~D"  p (lp-distance point1 point2 p)))))))
				
(defun test-position-of-best ()
  (trace position-of-best)
  (let ((seq '(1 2 3 4 5)))
    (dolist (aggr (list #'min #'max))
      (dolist (fitn (list (lambda (el) (- el))
			  (lambda (el) (expt el 2))
			  (lambda (el) (+ el 1))))
	(display 0 "~%....aggr=~A | fit~A~%~A~%" aggr fitn (position-of-best seq fitn aggr)))))
  (untrace position-of-best))  

(defun cluster-setup (max-val min-var var-range dist num-clusters num-points)
  "Sets globals used as defaults for generating cluster data"
  (setf *max-value* max-val
	*min-variance* min-var
	*variance-range* var-range
	*distribution* dist
	*num-clusters* num-clusters
	*num-points* num-points))
  
(defun generate-cluster (num-clusters num-points max-val min-var var-range dist)
  "Sets up cluster data, creates clusters, and prints them in csv format. Saves cluster data in *test* but returns nil."
  (cluster-setup max-val min-var var-range dist num-clusters num-points)
  (setf test (random-clusters))
  (csv-write-points test)
  nil)
