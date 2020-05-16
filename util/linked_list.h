/*     Chalmesh, 3-D overlapping grid generator */
/*     Copyright (C) 1997 Anders Petersson */

/*     This program is free software; you can redistribute it and/or modify */
/*     it under the terms of the GNU General Public License as published by */
/*     the Free Software Foundation; either version 2 of the License, or */
/*     (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/*     You can contact the author of Chalmesh by email at andersp@na.chalmers.se,  */
/*     or by paper-mail to  */

/*     Anders Petersson */
/*     Hydromechanics division */
/*     Department of Naval Architecture and Ocean Engineering */
/*     Chalmers University of Technology */
/*     412 96 Gothenburg */
/*     SWEDEN */
#ifndef linked_list_h
#define linked_list_h

typedef struct linked_list{
  struct linked_list_member *first, *last;
  int n_members;
} linked_list;

typedef struct linked_list_member {
  struct linked_list_member *prev, *next;
  void * data;
} linked_list_member;

linked_list *
new_linked_list(void);
linked_list *
delete_linked_list(linked_list *head);
linked_list_member *
new_link( linked_list *head );
linked_list_member *
new_link_before( linked_list *head, linked_list_member *old_link );
linked_list_member *
new_last_link( linked_list *head );
void 
delete_link( linked_list_member *this_link, linked_list *head );
void
link_insert_first( linked_list *head, linked_list_member *new_link );

#endif
