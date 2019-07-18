import { Component, OnInit } from '@angular/core';
import { NgForm } from '@angular/forms';
import { RequestService } from './requests.service';

@Component({
  selector: 'app-home',
  templateUrl: 'home.page.html',
  styleUrls: ['home.page.scss'],
})
export class HomePage implements OnInit {

  // Booleans that represents the current process. They have been used to show the relative ion-card
  public isTrimming: boolean;
  public isShortStack: boolean;
  public isMiRNA: boolean;
  public isPiRNA: boolean;
  public isOthers: boolean;
  public isNovelMirna: boolean;
  public isNovelPirna: boolean;

  private multiqc: any;


  constructor(private requestService: RequestService) {
    this.isTrimming = false;
    this.isShortStack = false;
    this.isPiRNA = false;
    this.isOthers = false;
    this.isMiRNA = false;
    this.isNovelMirna = false;
    this.isNovelPirna = false;
  }

  ngOnInit() {
    document.getElementById('card2').setAttribute('disabled', 'true');
    document.getElementById('card3').setAttribute('disabled', 'true');
    document.getElementById('card4').setAttribute('disabled', 'true');
    document.getElementById('card5').setAttribute('disabled', 'true');
    document.getElementById('card6').setAttribute('disabled', 'true');
    document.getElementById('card7').setAttribute('disabled', 'true');
  }

  private buildBody(form: NgForm) {
    return JSON.stringify({
      multimap: form.value['mismatches'],
      cores: form.value['core'],
      p_value: form.value['p_value'],
      p_value_adjusted: form.value['p_value_adjusted'],
      log_2_fold: form.value['log']
    });
  }

  onSubmit(form: NgForm, name: string) {
    if (form.valid) {
      if (name === 'f1') {
        const body = JSON.stringify({
          qual: form.value['qual'],
          adapter: form.value['adapter']
        });

        // this.multiqc = this.request.trimmingStep(body);
        document.getElementById('card2').setAttribute('disabled', 'false');
      } else if (name === 'f2') {
        const body = JSON.stringify({
          multimap: form.value['tresh'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });

        // this.requestService.shortStack(body)
        document.getElementById('card3').setAttribute('disabled', 'false');
      } else if (form.name === 'f3') {
        const body = this.buildBody(form);

        // this.requestService.mirnaDetection(body)
        document.getElementById('card4').setAttribute('disabled', 'false');
      } else if (name === 'f4') {
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card5').setAttribute('disabled', 'false');
      } else if (name === 'f5') {
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card6').setAttribute('disabled', 'false');
      } else if (name === 'f6') {
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card7').setAttribute('disabled', 'false');
      } else {
        const body = this.buildBody(form);
        // richiesta
      }
    }
  }

}
