
test_that("membership order", {
    mem <- c(4, 2, 4, 1, 3, 3, 1, 3)
    mem_new <- sort_membership(mem)
    expect_equal(as.integer(mem_new), c(3, 4, 3, 2, 1, 1, 2, 1))
    expect_equal((attr(mem_new, "order")), c(3, 1, 4, 2))

    mem <- c(4, 1, 1, 3, 2, 1, 2, 3)
    mem_new <- sort_membership(mem)
    expect_equal(as.integer(mem_new), mem)
    expect_equal((attr(mem_new, "order")), 1:4)
})
