/* A Linked-List Memory Sort
   by Philip J. Erdelsky
   pje@efgh.com
   http://www.alumni.caltech.edu/~pje/
*/

#include <stdio.h>
#include "../tdef.h"

void *sort_linked_list(void *p, unsigned index,
  int (*compare)(void *, void *, void *), void *pointer, unsigned long *pcount)
{
  unsigned base;
  unsigned long block_size;
	
	partptr g;	

  struct record
  {
    struct record *next[1];
    /* other members not directly accessed by this function */
  };

  struct tape
  {
    struct record *first, *last;
    unsigned long count;
  } tape[4];

  /* Distribute the records alternately to tape[0] and tape[1]. */

  tape[0].count = tape[1].count = 0L;
  tape[0].first = NULL;
  base = 0;
  while (p != NULL)
  {
    struct record *next = ((struct record *)p)->next[index];
    ((struct record *)p)->next[index] = tape[base].first;
    tape[base].first = ((struct record *)p);
    tape[base].count++;

    /*		g = *(partptr *)tape[base].first; */
    /*	fprintf(stderr,"tape[base].first = %g\n",g->pos[0]); */

		
    p = next;
    base ^= 1;
    /*fprintf(stderr,"base = %d (%d)\n",base,tape[base].count);*/
  }

  /* If the list is empty or contains only a single record, then */
  /* tape[1].count == 0L and this part is vacuous.               */

  for (base = 0, block_size = 1L; tape[base+1].count != 0L;
    base ^= 2, block_size <<= 1)
  {
    int dest;
    struct tape *tape0, *tape1;
    tape0 = tape + base;
    /*fprintf(stderr,"tape0 = %d\n",tape0);*/
    tape1 = tape + base + 1;
    /*fprintf(stderr,"tape1 = %d\n",tape1);*/
    dest = base ^ 2;
    /*fprintf(stderr,"dest = %d\n",dest);*/
    tape[dest].count = tape[dest+1].count = 0;
    for (; tape0->count != 0; dest ^= 1)
    {
      unsigned long n0, n1;
      struct tape *output_tape = tape + dest;
      n0 = n1 = block_size;
      while (1)
      {
        struct record *chosen_record;
        struct tape *chosen_tape;
        if (n0 == 0 || tape0->count == 0)
        {
          if (n1 == 0 || tape1->count == 0)
            break;
          chosen_tape = tape1;
          n1--;
        }
        else if (n1 == 0 || tape1->count == 0)
        {
          chosen_tape = tape0;
          n0--;
        }
        else if ((*compare)(tape0->first, tape1->first, pointer) > 0)
        {
          chosen_tape = tape1;
          n1--;
        }
        else
        {
          chosen_tape = tape0;
          n0--;
        }

	/* 				fprintf(stderr,"1 :: chosen_tape->count = %d\n",chosen_tape->count);*/
	      chosen_tape->count--;
	      /*				fprintf(stderr,"2 :: chosen_tape->count = %d\n",chosen_tape->count);*/
        chosen_record = chosen_tape->first;
        chosen_tape->first = chosen_record->next[index];
        if (output_tape->count == 0)
          output_tape->first = chosen_record;
        else
          output_tape->last->next[index] = chosen_record;
        output_tape->last = chosen_record;
        output_tape->count++;
      }
    }
  }

  if (tape[base].count > 1L)
    tape[base].last->next[index] = NULL;
  if (pcount != NULL)
    *pcount = tape[base].count;
  return tape[base].first;
}
