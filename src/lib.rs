use std::cmp::Ordering;
use std::convert::TryInto;
use std::mem;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    // use crate::coordinates_to_index;

    // #[test]
    // fn hilber_3() {
    //     let pt = coordinates_to_index(&vec![1, 2, 3], 8).unwrap();
    //     assert_eq!(pt, 36);
    // }
}

/// Represent an index on the Hilbert curve
pub type HilbertIndex = moore_hilbert_sys::BitmaskT;

/// Storage type for a single coordinate
pub type HilbertCoordinate = u64;

/// The datatype that contain the number of bits per dimension
pub type BitsPerDimensionType = usize;

/// The number of bits in each byte
const BITS_PER_BYTE: usize = 8;

/// Convert coordinates of a point on a Hilbert curve to its index.
/// # Arguments
///
/// * `bits_per_dimension` - Number of bits/coordinate.
/// * `coords` - Slice of coordinate values
///
/// # Returns
///
/// * `index` - Output index value.  nDims*nBits bits.
///
/// # Assumptions
///
/// `length of coords` * `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
/// ```
/// let bits_per_dimension = 8;
/// let r = moore_hilbert::coordinates_to_index(bits_per_dimension, &vec![1,2,3]).unwrap();
/// assert_eq!(r, 36);
/// ```
pub fn coordinates_to_index(
    bits_per_dimension: BitsPerDimensionType,
    coords: &[HilbertCoordinate],
) -> Result<HilbertIndex, ()> {
    if bits_per_dimension * coords.len() > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE {
        panic!("number of coordinates * bits_per_dimension > sizeof(HilbertIndex) * BITS_PER_BYTE");
    }

    unsafe {
        return Ok(moore_hilbert_sys::hilbert_c2i(
            coords.len().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            coords.as_ptr(),
        ));
    }
}

/// Convert an index into a Hilbert curve to a set of coordinates.
///
/// # Arguments
///
/// * `bits_per_dimension` - Number of bits per dimension
/// * `index` - The index, contains the number of dimensions * `nBits` bits (so `nDims` * `nBits` must be <= `BITS_PER_BYTE` * sizeof(`HilbertIndex`)).
/// * `coords` - The slice where the coordinates will be written
///
/// # Assumptions
///
/// `number of dimesions` * `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
/// # Example
///
/// ```
/// let bits_per_dimension = 8;
///
/// // Start by getting an index along the Hilbert curve for this point
/// let start_vec = vec![1,2,3];
/// let r = moore_hilbert::coordinates_to_index(bits_per_dimension, &start_vec).unwrap();
/// assert_eq!(r, 36);
///
/// // A place to put the coordinates from the Hilbert curve index
/// let mut extracted_coords = [0,0,0];
/// moore_hilbert::index_to_coordinates(bits_per_dimension, r, &mut extracted_coords);
///
/// /// The coordinates should match.
/// assert_eq!(start_vec, extracted_coords);
///
/// // increment the index and make sure the coords don't match
/// moore_hilbert::index_to_coordinates(bits_per_dimension, r+1, &mut extracted_coords);
/// assert_ne!(start_vec, extracted_coords);
///
/// ```
pub fn index_to_coordinates(
    bits_per_dimension: BitsPerDimensionType,
    index: HilbertIndex,
    coords: &mut [HilbertCoordinate],
) -> () {
    if bits_per_dimension as usize * coords.len() > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE {
        panic!("number of coordinates * bits_per_dimension > mem::size_of(HilbertIndex) * BITS_PER_BYTE");
    }

    unsafe {
        return moore_hilbert_sys::hilbert_i2c(
            coords.len().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            index,
            coords.as_mut_ptr(),
        );
    }
}

/// Determine which of two points lies further along the Hilbert curve
/// # Arguments
/// * `bits_per_dimension` - Number of bits/coordinate.
/// * `coord1` - Slice of nDims coordinates
/// * `coord2` - Slice of nDims coordinates
///
/// # Returns
///
/// * Ordering result that indicates the comparison of coord1 and coord2
///
/// # Assumptions
///
/// `nBits` <= (sizeof `BitmaskT`) * `bits_per_byte`
///
/// # Example
///
/// ```
/// use std::cmp::Ordering;
/// let coords = vec![
///   vec![1,2,3],
///   vec![1,2,4],
/// ];
///
/// let bits_per_dimension = 4;
///
/// assert_eq!(moore_hilbert::hilbert_cmp(bits_per_dimension, &coords[0], &coords[1]), Ordering::Less);
/// assert_eq!(moore_hilbert::hilbert_cmp(bits_per_dimension, &coords[0], &coords[0]), Ordering::Equal);
/// assert_eq!(moore_hilbert::hilbert_cmp(bits_per_dimension, &coords[1], &coords[0]), Ordering::Greater);
/// ```
pub fn hilbert_cmp(
    bits_per_dimension: BitsPerDimensionType,
    coord1: &[HilbertCoordinate],
    coord2: &[HilbertCoordinate],
) -> Ordering {
    if bits_per_dimension as usize > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE {
        panic!("bits_per_dimension > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE");
    }

    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    unsafe {
        let r = moore_hilbert_sys::hilbert_cmp(
            coord1.len().try_into().unwrap(),
            mem::size_of::<HilbertCoordinate>().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            coord1.as_ptr() as *const std::ffi::c_void,
            coord2.as_ptr() as *const std::ffi::c_void,
        );
        if r == -1 {
            return Ordering::Less;
        }
        if r == 0 {
            return Ordering::Equal;
        }
        return Ordering::Greater;
    }
}

/// Determine which of two points lies further along the Hilbert curve
///
/// # Arguments
///
/// * `coord1` - Array of nDims doubles
/// * `coord2` - Array of nDims doubles
///
/// # Returns
///
/// * Ordering result that indicates the comparison of coord1 and coord2
///
/// ```
/// use std::cmp::Ordering;
/// let coords = vec![
///   vec![1.0, 2.0, 3.0],
///   vec![1.0, 2.0, 4.0],
/// ];
///
///
/// assert_eq!(moore_hilbert::hilbert_ieee_cmp(&coords[0], &coords[0]), Ordering::Equal);
/// assert_eq!(moore_hilbert::hilbert_ieee_cmp(&coords[0], &coords[1]), Ordering::Greater);
/// assert_eq!(moore_hilbert::hilbert_ieee_cmp(&coords[1], &coords[0]), Ordering::Less);
/// ```
pub fn hilbert_ieee_cmp(coord1: &[f64], coord2: &[f64]) -> Ordering {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    unsafe {
        let r = moore_hilbert_sys::hilbert_ieee_cmp(
            coord1.len().try_into().unwrap(),
            coord1.as_ptr(),
            coord2.as_ptr(),
        );
        if r == -1 {
            return Ordering::Less;
        }
        if r == 0 {
            return Ordering::Equal;
        }
        return Ordering::Greater;
    }
}

/// Advance from one point to its successor on a Hilbert curve
///
/// # Arguments
///
/// * `bits_per_dimension` - Number of bits/coordinate.
/// * `coord` - Slice of coordinates changed to be the next point along the curve
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
/// ```
/// let bits_per_dimension = 8;
/// let mut coords = vec![2, 2, 5];
///
/// /// Get the initial position along the Hilbert curve
/// let first_index = moore_hilbert::coordinates_to_index(bits_per_dimension, &coords).unwrap();
///
/// /// Increment that position
/// moore_hilbert::hilbert_incr(bits_per_dimension, &mut coords);
///
/// /// Convert the incremented position back to a new index on the Hilbert curve
/// let new_index = moore_hilbert::coordinates_to_index(bits_per_dimension, &coords).unwrap();
///
/// /// The newly incremented index should advance along the curve by 1.
/// assert_eq!(new_index-first_index, 1);
///
///
/// ```
pub fn hilbert_incr(
    bits_per_dimension: BitsPerDimensionType,
    coord: &mut [HilbertCoordinate],
) -> () {
    if bits_per_dimension as usize > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE {
        panic!("bits_per_dimension > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE");
    }
    unsafe {
        moore_hilbert_sys::hilbert_incr(
            coord.len().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            coord.as_mut_ptr(),
        );
    }
}

/// Determine the first or last vertex of a box to lie on a Hilbert curve
///
/// # Arguments
/// * `bits_per_dimension`   - Number of bits per coordinate
/// * `find_min` - Is the least vertex sought?
/// * `c1`      - One corner of box
/// * `c2`      - Opposite corner
///
/// # Returns
///
/// `c1` and `c2` modified to refer to selected corner
/// value returned is log2 of size of largest power-of-two-aligned box that
/// contains the selected corner and no other corners
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
pub fn hilbert_box_vtx(
    bits_per_dimension: BitsPerDimensionType,
    find_min: bool,
    coord1: &mut [HilbertCoordinate],
    coord2: &mut [HilbertCoordinate],
) -> usize {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    if bits_per_dimension as usize > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE {
        panic!("bits_per_dimension > mem::size_of::<HilbertIndex>() * BITS_PER_BYTE");
    }

    unsafe {
        return moore_hilbert_sys::hilbert_box_vtx(
            coord1.len().try_into().unwrap(),
            mem::size_of::<HilbertCoordinate>().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            if find_min { 1 } else { 0 },
            coord1.as_ptr() as *mut std::ffi::c_void,
            coord2.as_ptr() as *mut std::ffi::c_void,
        ) as usize;
    }
}

/// Determine the first or last vertex of a box to lie on a Hilbert curve
///
/// # Arguments
///
/// * `find_min` - Is the least vertex sought?
/// * `c1`      - One corner of box
/// * `c2`      - Opposite corner
///
/// # Returns
/// `c1` and `c2` modified to refer to selected corder
/// value returned is log2 of size of largest power-of-two-aligned box that
/// contains the selected corner and no other corners
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
pub fn hilbert_ieee_box_vtx(find_min: bool, coord1: &mut [f64], coord2: &mut [f64]) -> usize {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    unsafe {
        return moore_hilbert_sys::hilbert_ieee_box_vtx(
            coord1.len().try_into().unwrap(),
            if find_min { 1 } else { 0 },
            coord1.as_mut_ptr(),
            coord2.as_mut_ptr(),
        ) as usize;
    }
}

/// Determine the first or last point of a box to lie on a Hilbert curve
///
/// # Arguments
///
/// * `bits_per_dimension` - Number of bits/coordinate.
/// * `find_min` - Is it the least vertex sought?
/// * `coord1` - Coordinates of one corner of box
/// * `coord2` - Coordinates of the opposite corner of box
///
/// # Returns
///
/// `coord1` and `coord2` are modified to refer to the least point
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
pub fn hilbert_box_pt(
    bits_per_dimension: BitsPerDimensionType,
    find_min: bool,
    coord1: &mut [HilbertCoordinate],
    coord2: &mut [HilbertCoordinate],
) -> usize {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    unsafe {
        return moore_hilbert_sys::hilbert_box_pt(
            coord1.len().try_into().unwrap(),
            mem::size_of::<HilbertCoordinate>().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            if find_min { 1 } else { 0 },
            coord1.as_mut_ptr() as *mut std::ffi::c_void,
            coord2.as_mut_ptr() as *mut std::ffi::c_void,
        ) as usize;
    }
}

/// Determine the first or last point of a box to lie on a Hilbert curve
///
/// # Arguments
///
/// * `bits_per_dimension` - Number of bits/coordinate.
/// * `find_min` - Is it the least vertex sought?
/// * `coord1` - Coordinates of one corner of box
/// * `coord2` - Coordinates of the opposite corner of box
///
/// # Returns
///
/// `coord1` and `coord2` are modified to refer to the least point
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex`) * (`bits_per_byte`)
///
pub fn hilbert_ieee_box_pt(find_min: bool, coord1: &mut [f64], coord2: &mut [f64]) -> usize {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    unsafe {
        return moore_hilbert_sys::hilbert_ieee_box_pt(
            coord1.len().try_into().unwrap(),
            if find_min { 1 } else { 0 },
            coord1.as_mut_ptr(),
            coord2.as_mut_ptr(),
        ) as usize;
    }
}

/// Determine the first point of a box after a given point to lie on a Hilbert curve
///
/// # Arguments
///
/// * `bits_per_dimension`    - Number of bits per dimension
/// * `find_prev`  -  Is the previous point sought?
/// * `coord1`    - Coordinates of  one corner of the box
/// * `coord2`    - Coordinates of the opposite corner of the box
/// * `point`    - Coordinates whish are a lower bound on point returned
///
/// # Returns
///
/// If true `coord1` and `coord2` are modified to point to the leeast point after point in the box
///
/// If false the arguments are unchanged and the point is beyond the last point of the box
///
/// # Assumptions
///
/// `bits_per_dimension` <= (sizeof `HilbertIndex) * (`bits_per_byte`)
///
pub fn hilbert_nextinbox(
    bits_per_dimension: BitsPerDimensionType,
    find_prev: bool,
    coord1: &mut [HilbertCoordinate],
    coord2: &mut [HilbertCoordinate],
    point: &[HilbertCoordinate],
) -> bool {
    if coord1.len() != coord2.len() {
        panic!("Coordinates supplied are not equal in length");
    }

    if point.len() != coord1.len() {
        panic!("Coordinates supplied are not equal in length to the point supplied");
    }

    unsafe {
        if moore_hilbert_sys::hilbert_nextinbox(
            coord1.len().try_into().unwrap(),
            mem::size_of::<HilbertCoordinate>().try_into().unwrap(),
            bits_per_dimension.try_into().unwrap(),
            if find_prev { 1 } else { 0 },
            coord1.as_mut_ptr() as *mut std::ffi::c_void,
            coord2.as_mut_ptr() as *mut std::ffi::c_void,
            point.as_ptr() as *mut std::ffi::c_void,
        ) == 1
        {
            true
        } else {
            false
        }
    }
}
