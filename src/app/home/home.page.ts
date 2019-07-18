import { Component } from '@angular/core';
import { NgForm } from '@angular/forms';

@Component({
  selector: 'app-home',
  templateUrl: 'home.page.html',
  styleUrls: ['home.page.scss'],
})
export class HomePage {

  constructor() { }

  onSubmit(form: NgForm) {
    if (form.valid) {
      if (form.name === 'f1') {
        // abilita la card n.2 dopo la richiesta
      } else if (form.name === 'f2') {
        // abilita la card n.3 dopo la richiesta
      } else if (form.name === 'f3') {
        // abilita la card n.4 dopo la richiesta
      } else if (form.name === 'f4') {
        // abilita la card n.5 dopo la richiesta
      } else if (form.name === 'f5') {
        // abilita la card n.6 dopo la richiesta
      } else if (form.name === 'f6') {
        // abilita la card n.7 dopo la richiesta
      } else { }

    }

  }

}
