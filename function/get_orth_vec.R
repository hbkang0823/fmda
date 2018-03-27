## Inputs ##
# obj: object for which points the orthogonal vectors will be found. P x 3.
# neighb_r: the radius in which the neighboring points are captured so that PCA can be performed.

## Outputs ##
# orthVec: orthogonal vectors.

get_orth_vec <- function(obj, neighb_r = 0.1) {

  r_pt <- colMeans(obj)	#reference point

  q_list <- 1:nrow(obj)

  q_vecmat <- matrix(NA, nrow=length(q_list), ncol=3)
  q_neighb_n <- rep(NA, length(q_list))

  for (q in q_list) {
    dist_to_q <- rep(0, 7150)
    for (p in 1:7150) {
      dist_to_q[p] <- sqrt(sum((obj[q,1:3]-obj[p,1:3])^2))
    }

    neighb_q <- obj[dist_to_q < neighb_r,]
    q_neighb_n[q] <- dim(neighb_q)[1]
    neighb_q.pca <- prcomp(neighb_q[,1:3])

    vec_q <- neighb_q.pca$rotation[,3]
    f_pt <- obj[q,1:3]		#facial point
    sign_vec <- sign(sum((f_pt - r_pt) * vec_q))	#determines sign
    q_vecmat[q,] <- sign_vec * vec_q

  }

  orthVec <- q_vecmat
  return(orthVec)
}
